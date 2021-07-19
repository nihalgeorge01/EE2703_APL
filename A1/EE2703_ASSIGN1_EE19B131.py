#! /usr/bin/env python3

# Author : Nihal John George (EE19B131)

'''
USAGE

$ python EE2703_ASSIGN1_EE19B131.py inp.txt
'''

import sys

CIRCUIT = '.circuit'
END = '.end'

# Base element class
class Element:
    def __init__(self, params):
        (self.name, 
        self.el_type, 
        self.n1, 
        self.n2, 
        self.n3, 
        self.n4, 
        self.v_cont, 
        self.value) = params
        
def token_parse(tokens):
    """
    Takes list of tokens and determines if they form a valid element, parses them into an Element object if valid

    Args:
        tokens : list(string)
    
    Returns:
        obj : Element class instance, or None if tokens are invalid
    """

    if len(tokens) < 4:
        return None
    else:
        name = tokens[0]
        el_type = tokens[0][0]
        n1 = tokens[1]
        n2 = tokens[2]
        n3 = None
        n4 = None
        v_cont = None
        value = None
        
        if not (n1.isalnum() and n2.isalnum()):                 # node names alphanumeric
            return None
        if tokens[0][0] in 'RLC' and len(tokens)==4:
            value = tokens[3]
        elif tokens[0][0] in 'VI' and len(tokens)==4:
            value = tokens[3]
        elif tokens[0][0] in 'EG' and len(tokens)==6:
            n3 = tokens[3]
            n4 = tokens[4]
            if not (n3.isalnum() and n4.isalnum()):             # node names alphanumeric
                return None
            value = tokens[5]
        elif tokens[0][0] in 'HF' and len(tokens)==5:
            v_cont = tokens[3]
            value = tokens[4]
        else:                                                   # Invalid line
            return None
        
        # Construct object with element data for future use
        params = [name, el_type, n1, n2, n3, n4, v_cont, value]
        return Element(params)

def netlist_to_data(path):
    """
    Takes path to netlist file and converts to data usable for circuit solving

    Args:
        path : string : path to file

    Returns:
        list of elements and their parameters
    """

    # Open file
    try:
        f = open(path, 'r')
    except FileNotFoundError:
        print("ERR: File not found")
        return None

    # Read file
    lines = f.readlines()
    f.close()

    # Remove comment portion, and extra leading and trailing space of each line
    lines = [i.split('#')[0].lstrip().rstrip() for i in lines]

    # Search for .circuit, exit on erroneous input
    start = None
    end = None
    for i in range(len(lines)):
        if lines[i][:len(CIRCUIT)] == CIRCUIT:
            if not start:
                start = i+1
            else:
                print("ERR: Extra .circuit found")
                return None
        elif lines[i][:len(END)] == END:
            if not start:
                print("ERR: .end found before .circuit")
                return None
            else:
                end = i
                break
    
    if not start:
        print("ERR: No .circuit found")
        return None
    if not end:
        print("ERR: No .end found")
        return None
    
    # Get all relevant lines, reverse the order of lines
    lines = lines[start:end][::-1]

    # Parse tokens into correct fields, exit on erroneous input
    obj_list= []
    for line in lines:
        
        # Omit empty lines
        if line == '':
            continue

        # Get space separated tokens
        tokens = line.split()
        
        # Assignment 1 Output
        print(' '.join(tokens[::-1]))

        # Parse tokens into objects
        el = token_parse(tokens)
        
        if not el:
            print("ERR: Invalid line here:", tokens)
            return None
    
        obj_list.append(el)
    
    return obj_list


if __name__ == '__main__':
    arg_count = len(sys.argv)
    if arg_count > 2:
        print("WARN: Too many arguments, considering only first one")
    elif arg_count < 2:
        print("ERR: Too less arguments")
        sys.exit(1)
    
    obj_list = netlist_to_data(sys.argv[1])

    # if None, some error has happened, so exit with code 1
    if obj_list == None:
        sys.exit(1)

    # normal exit code 0
    sys.exit(0)
    