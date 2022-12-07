
import numpy as np

# Accessory functions

# Name of input: value of input
# vars(args) -> returns dictionary from a namespace object
# so if args = parser.parse_args(), then vars(args) is the 
# dictionary of arguments : value pairs resulting from the parser

def dict_to_log(args_dict):
    """Create a formatted list of strings of the key value pairs of a dict."""
    args_list = [str(args_dict.keys[i])+
            " : {}".format(
                args_dict[args_dict.keys[i]]
                ) for i in range(len(args_dict))\
                if isinstance(args_dict[args_dict.keys[i]],
                (float, np.float32, np.float64))\
                else
                str(args_dict.keys[i])+
            " : {}".format(
                args_dict[args_dict.keys[i]]
                )]


def input_logger(filename, args_input):
    """Write argument input to file."""
    if len(filename) < 6 or filename[-6:] != ".input":
        filename = filename + ".input"
    
    if not isinstance(args_input, dict):
        args_input = vars(args_input)

    args_input_list = [str(args_input.keys[i])+
        " : {}".format(
            args_input[args_input.keys[i]]
            ) for i in range(len(args_input))]
    
    with open(filename, 'w') as f:
        f.writelines(args_input_list.sort())
    
    return


def data_logger(filename, data_dict):
    """Write (internal) data and parameters to file."""
    pass


def write_output_file(filename, args_input, data_dict):
    "Write output file with input args and data."
    input_logger(filename, args_input)
    data_logger(filename, data_dict)
    return


def dict_to_name(sep="_", **kwargs):
    """Construct a name from a dictionary"""
    name_as_list = list(map("".join, zip(kwargs.keys(), kwargs.values())))
    return sep.join(name_as_list)

def name_constructor(name="", sep="_", ):
    """Construct name dictionary, and return name"""
    pass
