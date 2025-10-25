from datetime import datetime

def cron(func, *args, **kwargs):
    """
    Decorator function to monitor the runtime of a function.
    
    This function takes another function as input, along with any number of positional and keyword arguments.
    It then defines a new function that wraps the input function, adding functionality to measure and print its runtime.
    
    Parameters:
    func (function): The function to be decorated.
    *args: Variable length argument list for the function to be decorated.
    **kwargs: Arbitrary keyword arguments for the function to be decorated.

    Returns:
    new_func (function): The decorated function with added functionality to print its runtime.
    """
    def new_func(*args, **kwargs):
        #print(f'Function "{func.__name__}" was called.')
        start_time = datetime.now()

        return_values = func(*args, **kwargs)

        end_time = datetime.now()
        run_time = end_time - start_time
        print(f'Runtime {func.__name__}: {run_time}\n')
        return return_values
    return new_func
