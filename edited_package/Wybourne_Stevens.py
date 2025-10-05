from .constants import WYBOURNE_STEVENS_CONSTS, LStheta, theta
from typing import Dict

# used outside the package

def WybourneToStevens(ion: str, Bdict: Dict[str, float], LS: bool = False) -> Dict[str, float]:
    """
    Convert Wybourne normalized parameters to Stevens notation.
    
    Wybourne and Stevens are two different normalization conventions for 
    crystal field parameters. This converts Wybourne B^k_q to Stevens B_n^m.
    
    Args:
        ion: Ion symbol (e.g., 'Yb3+', 'Dy3+')
        Bdict: Dictionary of Wybourne parameters with keys like 'B20', 'B40', etc.
            Values in meV.
        LS: If True, use LS coupling basis (for TM ions). If False, use J basis (RE ions).
    
    Returns:
        Dictionary of Stevens parameters with same keys, values in meV.
    
    Example:
        >>> wyb_params = {'B20': 0.5, 'B40': -0.02}
        >>> stev_params = WybourneToStevens('Yb3+', wyb_params)
        >>> print(stev_params['B20'])
    
    Note:
        Conversion: B_n^m (Stevens) = c_n^m * θ_n * B^n_m (Wybourne)
        where c_n^m from WYBOURNE_STEVENS_CONSTS, θ_n from theta() or LStheta()
    """
    return _convert_parameters(ion, Bdict, LS, to_stevens=True)


def StevensToWybourne(ion: str, Bdict: Dict[str, float], LS: bool = False) -> Dict[str, float]:
    """
    Convert Stevens parameters to Wybourne normalized notation.
    
    Inverse operation of WybourneToStevens. Converts Stevens B_n^m to Wybourne B^k_q.
    
    Args:
        ion: Ion symbol (e.g., 'Yb3+', 'Dy3+')
        Bdict: Dictionary of Stevens parameters with keys like 'B20', 'B40', etc.
            Values in meV.
        LS: If True, use LS coupling basis (for TM ions). If False, use J basis (RE ions).
    
    Returns:
        Dictionary of Wybourne parameters with same keys, values in meV.
    
    Example:
        >>> stev_params = {'B20': 0.5, 'B40': -0.02}
        >>> wyb_params = StevensToWybourne('Yb3+', stev_params)
        >>> print(wyb_params['B20'])
    
    Note:
        Conversion: B^n_m (Wybourne) = B_n^m (Stevens) / (c_n^m * θ_n)
    """
    return _convert_parameters(ion, Bdict, LS, to_stevens=False)


def _convert_parameters(ion: str, Bdict: Dict[str, float], LS: bool, 
                       to_stevens: bool) -> Dict[str, float]:
    """
    Helper function to convert between Wybourne and Stevens parameters.
    
    Args:
        ion: Ion symbol
        Bdict: Input parameter dictionary
        LS: Whether to use LS coupling (True) or J basis (False)
        to_stevens: If True, convert Wyb→Stevens. If False, convert Stevens→Wyb.
    
    Returns:
        Converted parameter dictionary
    """
    result = {}
    theta_func = LStheta if LS else theta
    
    for param_key in Bdict:             # 'B42' → n=4, m=2;
        n = int(param_key[1])
        m = int(param_key[2:])
        
        conversion_factor = WYBOURNE_STEVENS_CONSTS[n][m] * theta_func(ion, n)
        
        if to_stevens:
            result['B' + param_key[1:]] = Bdict[param_key] * conversion_factor
        else:
            result['B' + param_key[1:]] = Bdict[param_key] / conversion_factor
    
    return result