import ast
import os
import sys
from pathlib import Path
from collections import defaultdict, namedtuple
from typing import Dict, List, Set, Tuple, Optional, Any
import json

# Data structures for tracking entities and their usage
Entity = namedtuple('Entity', ['type', 'name', 'file', 'line', 'parent'])
Usage = namedtuple('Usage', ['entity', 'used_in_file', 'context', 'line'])

class EntityVisitor(ast.NodeVisitor):
    """AST visitor to extract entities from a Python file."""
    
    def __init__(self, filename: str):
        self.filename = filename
        self.entities = []
        self.current_class = None
        self.current_function = None
        self.imports = {}  # Track imports for resolution
        
    def visit_ClassDef(self, node: ast.ClassDef):
        """Extract class definitions."""
        self.entities.append(Entity(
            type='CLASS',
            name=node.name,
            file=self.filename,
            line=node.lineno,
            parent=None
        ))
        
        # Track current class for method context
        old_class = self.current_class
        self.current_class = node.name
        
        # Visit methods within the class
        for item in node.body:
            if isinstance(item, ast.FunctionDef):
                self.entities.append(Entity(
                    type='METHOD',
                    name=item.name,
                    file=self.filename,
                    line=item.lineno,
                    parent=node.name
                ))
        
        self.generic_visit(node)
        self.current_class = old_class
        
    def visit_FunctionDef(self, node: ast.FunctionDef):
        """Extract function definitions."""
        if self.current_class is None:  # Only standalone functions
            self.entities.append(Entity(
                type='FUNCTION',
                name=node.name,
                file=self.filename,
                line=node.lineno,
                parent=None
            ))
            
            old_function = self.current_function
            self.current_function = node.name
            self.generic_visit(node)
            self.current_function = old_function
        else:
            self.generic_visit(node)
    
    def visit_Assign(self, node: ast.Assign):
        """Extract module-level variables and constants."""
        if self.current_class is None and self.current_function is None:
            for target in node.targets:
                if isinstance(target, ast.Name):
                    # Check if it's a constant (UPPER_CASE convention)
                    entity_type = 'CONSTANT' if target.id.isupper() else 'VARIABLE'
                    self.entities.append(Entity(
                        type=entity_type,
                        name=target.id,
                        file=self.filename,
                        line=node.lineno,
                        parent=None
                    ))
        self.generic_visit(node)
        
    def visit_Import(self, node: ast.Import):
        """Track import statements."""
        for alias in node.names:
            self.imports[alias.asname or alias.name] = alias.name
        self.generic_visit(node)
        
    def visit_ImportFrom(self, node: ast.ImportFrom):
        """Track from-import statements."""
        module = node.module or ''
        for alias in node.names:
            imported_name = alias.name
            local_name = alias.asname or imported_name
            self.imports[local_name] = f"{module}.{imported_name}" if module else imported_name
        self.generic_visit(node)


class UsageVisitor(ast.NodeVisitor):
    """AST visitor to find entity usage within a file."""
    
    def __init__(self, filename: str, all_entities: Dict[str, List[Entity]]):
        self.filename = filename
        self.all_entities = all_entities
        self.usages = []
        self.current_class = None
        self.current_method = None
        self.current_function = None
        self.imports = {}
        self.imported_modules = set()
        
    def get_context(self) -> str:
        """Get the current context string."""
        if self.current_class and self.current_method:
            return f"class:{self.current_class} → method:{self.current_method}"
        elif self.current_class:
            return f"class:{self.current_class}"
        elif self.current_function:
            return f"function:{self.current_function}"
        else:
            return "module-level"
    
    def visit_ClassDef(self, node: ast.ClassDef):
        """Track class context."""
        old_class = self.current_class
        self.current_class = node.name
        self.generic_visit(node)
        self.current_class = old_class
        
    def visit_FunctionDef(self, node: ast.FunctionDef):
        """Track function/method context."""
        if self.current_class:
            old_method = self.current_method
            self.current_method = node.name
            self.generic_visit(node)
            self.current_method = old_method
        else:
            old_function = self.current_function
            self.current_function = node.name
            self.generic_visit(node)
            self.current_function = old_function
            
    def visit_Name(self, node: ast.Name):
        """Track usage of names (variables, functions, classes)."""
        name = node.id
        
        # Check if this name matches any entity
        for entity_list in self.all_entities.values():
            for entity in entity_list:
                if entity.name == name and entity.file != self.filename:
                    # Found usage of an entity from another file
                    self.usages.append(Usage(
                        entity=entity,
                        used_in_file=self.filename,
                        context=self.get_context(),
                        line=node.lineno
                    ))
        
        self.generic_visit(node)
        
    def visit_Attribute(self, node: ast.Attribute):
        """Track attribute access (e.g., module.function)."""
        # Handle chained attributes
        if isinstance(node.value, ast.Name):
            base_name = node.value.id
            attr_name = node.attr
            
            # Check if this matches any entity
            for entity_list in self.all_entities.values():
                for entity in entity_list:
                    if entity.name == attr_name:
                        self.usages.append(Usage(
                            entity=entity,
                            used_in_file=self.filename,
                            context=self.get_context(),
                            line=node.lineno
                        ))
        
        self.generic_visit(node)
        
    def visit_Call(self, node: ast.Call):
        """Track function/method calls."""
        # Handle direct function calls
        if isinstance(node.func, ast.Name):
            func_name = node.func.id
            
            for entity_list in self.all_entities.values():
                for entity in entity_list:
                    if entity.name == func_name and entity.type in ('FUNCTION', 'CLASS', 'METHOD'):
                        if entity.file != self.filename:
                            self.usages.append(Usage(
                                entity=entity,
                                used_in_file=self.filename,
                                context=self.get_context(),
                                line=node.lineno
                            ))
        
        self.generic_visit(node)
        
    def visit_Import(self, node: ast.Import):
        """Track imports."""
        for alias in node.names:
            self.imports[alias.asname or alias.name] = alias.name
            self.imported_modules.add(alias.name)
        self.generic_visit(node)
        
    def visit_ImportFrom(self, node: ast.ImportFrom):
        """Track from imports."""
        module = node.module or ''
        for alias in node.names:
            imported_name = alias.name
            local_name = alias.asname or imported_name
            self.imports[local_name] = imported_name
            
            # Check if this import matches any entity
            for entity_list in self.all_entities.values():
                for entity in entity_list:
                    if entity.name == imported_name:
                        self.usages.append(Usage(
                            entity=entity,
                            used_in_file=self.filename,
                            context="import-statement",
                            line=node.lineno
                        ))
        
        self.generic_visit(node)


class PythonPackageAnalyzer:
    """Main analyzer for Python package dependencies and usage."""
    
    def __init__(self, package_path: str):
        self.package_path = Path(package_path)
        self.entities: Dict[str, List[Entity]] = defaultdict(list)
        self.usages: List[Usage] = []
        self.file_contents: Dict[str, str] = {}
        
    def analyze(self):
        """Perform complete analysis of the package."""
        print(f"🔍 Analyzing package: {self.package_path}")
        print("=" * 80)
        
        # Phase 1: Discover all entities
        self._discover_entities()
        
        # Phase 2: Track usage across files
        self._track_usage()
        
        # Phase 3: Generate report
        self._generate_report()
        
    def _discover_entities(self):
        """Discover all entities in Python files."""
        python_files = list(self.package_path.rglob("*.py"))
        
        print(f"\n📂 Found {len(python_files)} Python files")
        print("-" * 40)
        
        for py_file in python_files:
            if "__pycache__" in str(py_file):
                continue
                
            try:
                with open(py_file, 'r', encoding='utf-8') as f:
                    content = f.read()
                    self.file_contents[str(py_file)] = content
                    
                tree = ast.parse(content, filename=str(py_file))
                visitor = EntityVisitor(str(py_file))
                visitor.visit(tree)
                
                self.entities[str(py_file)] = visitor.entities
                
                if visitor.entities:
                    print(f"  ✓ {py_file.name}: {len(visitor.entities)} entities found")
                    
            except Exception as e:
                print(f"  ✗ Error parsing {py_file.name}: {e}")
                
    def _track_usage(self):
        """Track usage of entities across all files."""
        print(f"\n🔗 Tracking cross-file dependencies...")
        print("-" * 40)
        
        for py_file, content in self.file_contents.items():
            try:
                tree = ast.parse(content, filename=py_file)
                visitor = UsageVisitor(py_file, self.entities)
                visitor.visit(tree)
                
                self.usages.extend(visitor.usages)
                
                if visitor.usages:
                    print(f"  ✓ {Path(py_file).name}: {len(visitor.usages)} dependencies found")
                    
            except Exception as e:
                print(f"  ✗ Error analyzing {Path(py_file).name}: {e}")
                
    def _generate_report(self):
        """Generate comprehensive usage report."""
        print("\n" + "=" * 80)
        print("📊 COMPREHENSIVE DEPENDENCY ANALYSIS REPORT")
        print("=" * 80)
        
        # Group usages by entity
        usage_by_entity: Dict[Entity, List[Usage]] = defaultdict(list)
        for usage in self.usages:
            usage_by_entity[usage.entity].append(usage)
        
        # Report by file
        for file_path, entities in self.entities.items():
            if not entities:
                continue
                
            rel_path = Path(file_path).relative_to(self.package_path)
            print(f"\n📁 {rel_path}")
            print("─" * 60)
            
            # Group entities by type
            by_type = defaultdict(list)
            for entity in entities:
                by_type[entity.type].append(entity)
            
            # Display classes and their methods
            for class_entity in by_type.get('CLASS', []):
                usages = usage_by_entity.get(class_entity, [])
                usage_count = len(usages)
                
                print(f"\n  🏛️  CLASS: {class_entity.name} (line {class_entity.line}) " +
                      f"[used by {usage_count} {'file' if usage_count == 1 else 'files'}]")
                
                # Show methods of this class
                class_methods = [e for e in by_type.get('METHOD', []) 
                                if e.parent == class_entity.name]
                for method in class_methods:
                    method_usages = usage_by_entity.get(method, [])
                    print(f"      └── {method.name}() (line {method.line}) " +
                          f"[used {len(method_usages)} times]")
                    
                    # Show usage details for methods
                    if method_usages:
                        for usage in method_usages[:3]:  # Show first 3 usages
                            rel_usage_path = Path(usage.used_in_file).relative_to(self.package_path)
                            print(f"          └── 📦 {rel_usage_path}:")
                            print(f"              └── {usage.context} (line {usage.line})")
                        if len(method_usages) > 3:
                            print(f"          └── ... and {len(method_usages) - 3} more usages")
                
                # Show usage details for class
                if usages:
                    print(f"      📍 Class usage locations:")
                    for usage in usages[:3]:
                        rel_usage_path = Path(usage.used_in_file).relative_to(self.package_path)
                        print(f"          └── 📦 {rel_usage_path}:")
                        print(f"              └── {usage.context} (line {usage.line})")
                    if len(usages) > 3:
                        print(f"          └── ... and {len(usages) - 3} more usages")
            
            # Display standalone functions
            for func_entity in by_type.get('FUNCTION', []):
                usages = usage_by_entity.get(func_entity, [])
                usage_count = len(usages)
                
                print(f"\n  📊 FUNCTION: {func_entity.name}() (line {func_entity.line}) " +
                      f"[used by {usage_count} {'file' if usage_count == 1 else 'files'}]")
                
                if usages:
                    for usage in usages[:3]:
                        rel_usage_path = Path(usage.used_in_file).relative_to(self.package_path)
                        print(f"      └── 📦 {rel_usage_path}:")
                        print(f"          └── {usage.context} (line {usage.line})")
                    if len(usages) > 3:
                        print(f"      └── ... and {len(usages) - 3} more usages")
            
            # Display constants
            for const_entity in by_type.get('CONSTANT', []):
                usages = usage_by_entity.get(const_entity, [])
                usage_count = len(usages)
                
                print(f"\n  🔧 CONSTANT: {const_entity.name} (line {const_entity.line}) " +
                      f"[used by {usage_count} {'file' if usage_count == 1 else 'files'}]")
                
                if usages:
                    for usage in usages[:2]:
                        rel_usage_path = Path(usage.used_in_file).relative_to(self.package_path)
                        print(f"      └── 📦 {rel_usage_path}:")
                        print(f"          └── {usage.context} (line {usage.line})")
                    if len(usages) > 2:
                        print(f"      └── ... and {len(usages) - 2} more usages")
            
            # Display variables
            for var_entity in by_type.get('VARIABLE', []):
                usages = usage_by_entity.get(var_entity, [])
                if usages:  # Only show variables that are actually used
                    print(f"\n  📌 VARIABLE: {var_entity.name} (line {var_entity.line}) " +
                          f"[used by {len(usages)} {'file' if len(usages) == 1 else 'files'}]")
        
        # Summary statistics
        self._print_summary_statistics(usage_by_entity)
        
    def _print_summary_statistics(self, usage_by_entity: Dict[Entity, List[Usage]]):
        """Print summary statistics about the package."""
        print("\n" + "=" * 80)
        print("📈 SUMMARY STATISTICS")
        print("=" * 80)
        
        # Total counts
        total_entities = sum(len(entities) for entities in self.entities.values())
        total_files = len(self.entities)
        total_usages = len(self.usages)
        
        print(f"\n  📊 Total files analyzed: {total_files}")
        print(f"  📊 Total entities discovered: {total_entities}")
        print(f"  📊 Total cross-file dependencies: {total_usages}")
        
        # Most used entities
        most_used = sorted(usage_by_entity.items(), key=lambda x: len(x[1]), reverse=True)[:5]
        if most_used:
            print("\n  🏆 Most frequently used entities:")
            for entity, usages in most_used:
                rel_path = Path(entity.file).relative_to(self.package_path)
                print(f"      {entity.type}: {entity.name} ({rel_path}) - {len(usages)} usages")
        
        # Files with most dependencies
        file_dependency_count = defaultdict(int)
        for usage in self.usages:
            file_dependency_count[usage.used_in_file] += 1
        
        most_dependent = sorted(file_dependency_count.items(), key=lambda x: x[1], reverse=True)[:3]
        if most_dependent:
            print("\n  🔗 Files with most dependencies:")
            for file_path, count in most_dependent:
                rel_path = Path(file_path).relative_to(self.package_path)
                print(f"      {rel_path}: {count} dependencies")
        
        # Unused entities (potentially dead code)
        unused_entities = []
        for entities in self.entities.values():
            for entity in entities:
                if entity not in usage_by_entity and entity.type != 'METHOD':
                    # Exclude methods as they might be called dynamically
                    unused_entities.append(entity)
        
        if unused_entities:
            print(f"\n  ⚠️  Potentially unused entities: {len(unused_entities)}")
            for entity in unused_entities[:5]:
                rel_path = Path(entity.file).relative_to(self.package_path)
                print(f"      {entity.type}: {entity.name} ({rel_path})")
            if len(unused_entities) > 5:
                print(f"      ... and {len(unused_entities) - 5} more")
    
    def export_to_json(self, output_file: str = "dependency_analysis.json"):
        """Export analysis results to JSON format."""
        export_data = {
            "package_path": str(self.package_path),
            "entities": [],
            "usages": [],
            "statistics": {}
        }
        
        # Convert entities to dictionaries
        for file_path, entities in self.entities.items():
            for entity in entities:
                export_data["entities"].append({
                    "type": entity.type,
                    "name": entity.name,
                    "file": file_path,
                    "line": entity.line,
                    "parent": entity.parent
                })
        
        # Convert usages to dictionaries
        for usage in self.usages:
            export_data["usages"].append({
                "entity_name": usage.entity.name,
                "entity_type": usage.entity.type,
                "entity_file": usage.entity.file,
                "used_in_file": usage.used_in_file,
                "context": usage.context,
                "line": usage.line
            })
        
        # Add statistics
        export_data["statistics"] = {
            "total_files": len(self.entities),
            "total_entities": sum(len(entities) for entities in self.entities.values()),
            "total_dependencies": len(self.usages)
        }
        
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(export_data, f, indent=2)
        
        print(f"\n💾 Analysis exported to: {output_file}")
    
    def export_to_markdown(self, output_file: str = "dependency_analysis.md"):
        """Export analysis results to Markdown with collapsible sections."""
        from datetime import datetime
        
        # Group usages by entity
        usage_by_entity: Dict[Entity, List[Usage]] = defaultdict(list)
        for usage in self.usages:
            usage_by_entity[usage.entity].append(usage)
        
        md_lines = []
        md_lines.append(f"# 📊 Python Package Dependency Analysis Report")
        md_lines.append(f"\n**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        md_lines.append(f"\n**Package:** `{self.package_path}`\n")
        md_lines.append("---\n")
        
        # Generate table of contents
        md_lines.append("## 📑 Table of Contents\n")
        for file_path in self.entities.keys():
            rel_path = Path(file_path).relative_to(self.package_path)
            anchor = str(rel_path).replace('/', '-').replace('\\', '-').replace('.', '-')
            md_lines.append(f"- [{rel_path}](#{anchor})")
        md_lines.append("\n---\n")
        
        # Report by file
        for file_path, entities in self.entities.items():
            if not entities:
                continue
                
            rel_path = Path(file_path).relative_to(self.package_path)
            anchor = str(rel_path).replace('/', '-').replace('\\', '-').replace('.', '-')
            md_lines.append(f"## 📁 {rel_path} {{#{anchor}}}\n")
            
            # Group entities by type
            by_type = defaultdict(list)
            for entity in entities:
                by_type[entity.type].append(entity)
            
            # Classes and methods
            if by_type.get('CLASS'):
                md_lines.append("### 🏛️ Classes\n")
                for class_entity in by_type['CLASS']:
                    usages = usage_by_entity.get(class_entity, [])
                    md_lines.append(f"#### **{class_entity.name}** *(line {class_entity.line})*\n")
                    md_lines.append(f"- **Used by:** {len(usages)} {'file' if len(usages) <= 1 else 'files'}\n")
                    
                    # Methods
                    class_methods = [e for e in by_type.get('METHOD', []) 
                                if e.parent == class_entity.name]
                    if class_methods:
                        md_lines.append("\n**Methods:**\n")
                        for method in class_methods:
                            method_usages = usage_by_entity.get(method, [])
                            md_lines.append(f"- `{method.name}()` *(line {method.line})* - {len(method_usages)} usages")
                            
                            if method_usages:
                                if len(method_usages) <= 3:
                                    md_lines.append("  <details><summary>View usages</summary>\n")
                                else:
                                    md_lines.append(f"  <details><summary>View all {len(method_usages)} usages</summary>\n")
                                
                                for usage in method_usages:
                                    rel_usage_path = Path(usage.used_in_file).relative_to(self.package_path)
                                    md_lines.append(f"  - `{rel_usage_path}`: {usage.context} (line {usage.line})")
                                md_lines.append("  </details>\n")
                    
                    # Class usages
                    if usages:
                        md_lines.append(f"\n<details><summary>📍 View all {len(usages)} class usages</summary>\n")
                        for usage in usages:
                            rel_usage_path = Path(usage.used_in_file).relative_to(self.package_path)
                            md_lines.append(f"- `{rel_usage_path}`: {usage.context} (line {usage.line})")
                        md_lines.append("</details>\n")
            
            # Standalone functions
            if by_type.get('FUNCTION'):
                md_lines.append("### 📊 Functions\n")
                for func_entity in by_type['FUNCTION']:
                    usages = usage_by_entity.get(func_entity, [])
                    md_lines.append(f"#### **{func_entity.name}()** *(line {func_entity.line})*\n")
                    md_lines.append(f"- **Used by:** {len(usages)} {'file' if len(usages) <= 1 else 'files'}\n")
                    
                    if usages:
                        md_lines.append(f"<details><summary>View all {len(usages)} usages</summary>\n")
                        for usage in usages:
                            rel_usage_path = Path(usage.used_in_file).relative_to(self.package_path)
                            md_lines.append(f"- `{rel_usage_path}`: {usage.context} (line {usage.line})")
                        md_lines.append("</details>\n")
            
            # Constants
            if by_type.get('CONSTANT'):
                md_lines.append("### 🔧 Constants\n")
                for const_entity in by_type['CONSTANT']:
                    usages = usage_by_entity.get(const_entity, [])
                    md_lines.append(f"- **{const_entity.name}** *(line {const_entity.line})* - {len(usages)} usages")
                    
                    if usages:
                        md_lines.append(f"  <details><summary>View usages</summary>\n")
                        for usage in usages:
                            rel_usage_path = Path(usage.used_in_file).relative_to(self.package_path)
                            md_lines.append(f"  - `{rel_usage_path}`: {usage.context} (line {usage.line})")
                        md_lines.append("  </details>\n")
            
            # Variables (only if used)
            used_vars = [v for v in by_type.get('VARIABLE', []) if usage_by_entity.get(v)]
            if used_vars:
                md_lines.append("### 📌 Variables\n")
                for var_entity in used_vars:
                    usages = usage_by_entity.get(var_entity, [])
                    md_lines.append(f"- **{var_entity.name}** *(line {var_entity.line})* - {len(usages)} usages")
                    
                    if usages:
                        md_lines.append(f"  <details><summary>View usages</summary>\n")
                        for usage in usages:
                            rel_usage_path = Path(usage.used_in_file).relative_to(self.package_path)
                            md_lines.append(f"  - `{rel_usage_path}`: {usage.context} (line {usage.line})")
                        md_lines.append("  </details>\n")
            
            md_lines.append("---\n")
        
        # Summary statistics
        md_lines.append("## 📈 Summary Statistics\n")
        
        total_entities = sum(len(entities) for entities in self.entities.values())
        total_files = len(self.entities)
        total_usages = len(self.usages)
        
        md_lines.append(f"- **Total files analyzed:** {total_files}")
        md_lines.append(f"- **Total entities discovered:** {total_entities}")
        md_lines.append(f"- **Total cross-file dependencies:** {total_usages}\n")
        
        # Most used entities
        most_used = sorted(usage_by_entity.items(), key=lambda x: len(x[1]), reverse=True)[:10]
        if most_used:
            md_lines.append("### 🏆 Most Frequently Used Entities\n")
            md_lines.append("| Entity | Type | File | Usage Count |")
            md_lines.append("|--------|------|------|-------------|")
            for entity, usages in most_used:
                rel_path = Path(entity.file).relative_to(self.package_path)
                md_lines.append(f"| `{entity.name}` | {entity.type} | {rel_path} | {len(usages)} |")
        
        # Files with most dependencies
        file_dependency_count = defaultdict(int)
        for usage in self.usages:
            file_dependency_count[usage.used_in_file] += 1
        
        most_dependent = sorted(file_dependency_count.items(), key=lambda x: x[1], reverse=True)[:5]
        if most_dependent:
            md_lines.append("\n### 🔗 Files with Most Dependencies\n")
            md_lines.append("| File | Dependency Count |")
            md_lines.append("|------|------------------|")
            for file_path, count in most_dependent:
                rel_path = Path(file_path).relative_to(self.package_path)
                md_lines.append(f"| {rel_path} | {count} |")
        
        # Unused entities
        unused_entities = []
        for entities in self.entities.values():
            for entity in entities:
                if entity not in usage_by_entity and entity.type != 'METHOD':
                    unused_entities.append(entity)
        
        if unused_entities:
            md_lines.append(f"\n### ⚠️ Potentially Unused Entities ({len(unused_entities)} total)\n")
            md_lines.append("<details><summary>View all unused entities</summary>\n")
            md_lines.append("| Entity | Type | File | Line |")
            md_lines.append("|--------|------|------|------|")
            for entity in unused_entities:
                rel_path = Path(entity.file).relative_to(self.package_path)
                md_lines.append(f"| `{entity.name}` | {entity.type} | {rel_path} | {entity.line} |")
            md_lines.append("\n</details>")
        
        # Write to file
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write('\n'.join(md_lines))
        
        print(f"\n📝 Markdown report exported to: {output_file}")



def main():
    """Main entry point."""
    if len(sys.argv) < 2:
        print("Usage: python analyzer.py <package_path> [--export] [--markdown]")
        print("\nExample: python analyzer.py ./my_package --markdown")
        sys.exit(1)
    
    package_path = sys.argv[1]
    
    if not os.path.exists(package_path):
        print(f"Error: Path '{package_path}' does not exist")
        sys.exit(1)
    
    if not os.path.isdir(package_path):
        print(f"Error: Path '{package_path}' is not a directory")
        sys.exit(1)
    
    # Run analysis
    analyzer = PythonPackageAnalyzer(package_path)
    analyzer.analyze()
    
    # Check for export options
    if "--export" in sys.argv:
        analyzer.export_to_json()
    
    if "--markdown" in sys.argv:
        analyzer.export_to_markdown()


if __name__ == "__main__":
    main()