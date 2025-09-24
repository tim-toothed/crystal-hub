# Отредактированный пакет PyCrystalFields

Используй данную команду в терминале, чтобы получить markdown с подробным анализом dependencies (где какие объекты прописаны и используются): 

```powershell
python analyzer.py edited_package --markdown
```

Вывод:
```markdown
📁 stevens_operators.py
────────────────────────────────────────────────────────────

  📊 FUNCTION: StevensOp() (line 5) [used by 14 files]
      └── 📦 cf_levels.py:
          └── import-statement (line 11)
      └── 📦 cf_levels.py:
          └── class:CFLevels → method:Bdict (line 47)
      └── 📦 cf_levels.py:
          └── class:CFLevels → method:Bdict (line 47)
      └── ... and 11 more usages

  📊 FUNCTION: LS_StevensOp() (line 102) [used by 10 files]
      └── 📦 cf_levels.py:
          └── import-statement (line 11)
      └── 📦 cf_levels.py:
          └── class:LS_CFLevels → method:Bdict (line 840)
      └── 📦 cf_levels.py:
          └── class:LS_CFLevels → method:Bdict (line 840)
      └── ... and 7 more usages
```

Результат - [dependency_analysis.md](../dependency_analysis.md)