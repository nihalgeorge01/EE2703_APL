#! /usr/bin/env python3

# Author : Nihal John George (EE19B131)

'''
USAGE

$ python EE2703_ASSIGN2_EE19B131.py ckt.netlist
'''

import sys
import numpy as np
ffs = np.format_float_scientific

CIRCUIT = '.circuit'
END = '.end'
AC = '.ac'
EPS = 1e-30                     # Infinitesimal
INF = 1/EPS                     # Infinity

# Base element class
class Element:
    def __init__(self, params):
        (self.name, 
        self.el_type,           # R L C V I E G H F
        self.n1, 
        self.n2, 
        self.form,              # ac, dc
        self.freq,              # frequency in Hz (impedance calculated using w=2*pi*freq)
        self.n3, 
        self.n4, 
        self.v_cont,            # controlling current passing through voltage source
        self.value,             # component value
        self.imp,               # complex impedance
        self.phase) = params
        
def token_parse(tokens):
    """
    Takes list of tokens and determines if they form a valid element, 
    parses them into an Element object if valid

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
        form = None
        freq = None
        n3 = None
        n4 = None
        v_cont = None
        value = None
        imp = None
        phase = None
        
        if not (n1.isalnum() and n2.isalnum()):                 # node names alphanumeric
            print("ERR: Node names must be alphanumeric:", n1, n2)
            return None

        if n1 == n2:
            print("ERR: Connecting nodes must be distinct: ", tokens)
            return None
        
        if tokens[0][0] in 'RLC' and len(tokens)==4:
            value = tokens[3]
            if tokens[0][0] == 'R':
                imp = value

            if n1 == 'GND':
                n1,n2 = n2,n1

        elif tokens[0][0] in 'VI' and len(tokens)==5 and tokens[3]=='dc':
            form = 'dc'
            value = tokens[4]

        elif tokens[0][0] in 'VI' and len(tokens)==6 and tokens[3]=='ac':
            form = 'ac'
            value = tokens[4]
            phase = tokens[5]

        elif tokens[0][0] in 'EG' and len(tokens)==6:
            n3 = tokens[3]
            n4 = tokens[4]
            if not (n3.isalnum() and n4.isalnum()):             # node names alphanumeric
                print("ERR: Node names must be alphanumeric:", n3, n4)
                return None
            if n3 == n4:
                print("ERR: Voltage-controlled source received same nodes for control:", n3, n4)
            value = tokens[5]

        elif tokens[0][0] in 'HF' and len(tokens)==5:
            v_cont = tokens[3]
            value = tokens[4]

        else:                                                   # Invalid line
            print("ERR: Invalid line tokens:", tokens)
            return None
        
        quants = [freq, value, phase, imp]
        for ind in range(len(quants)):
            q = quants[ind]
            e_q = None

            if q != None:
                e_q = get_quant(q)
                if e_q == None:
                    print("ERR: Invalid quantity: ", q)
                    return None
            
            quants[ind] = e_q

        freq, value, phase, imp = quants

        if tokens[3] == 'ac':
            value /= 2.0                                        # amplitude = Vpp/2
        
        # Construct object with element data for future use
        params = [name, el_type, n1, n2, form, freq, n3, n4, v_cont, value, imp, phase]
        return Element(params)

def get_quant(q):
    """
    Takes a string and returns it in a numeric type if suitable (using eval())

    Args:
        q : string : String to convert into numeric quantity

    Returns:
        e_q : numeric : Numeric quantity, or None if not suitable string
    """

    try:
        e_q = eval(q)
    except:
        return None

    if isinstance(e_q, (int,float,complex)):
        return e_q
    
    return None
    
def ac_parse(obj_list, ac_lines):
    """
    Parses the .ac lines and modifies the frequency of each voltage source
    and set impedances of passive elements
    WARN: Currently supports only single frequency
    
    Args:
        obj_list : list(Element) : list of elements
        ac_lines : list(string) : list of .ac commands

    Returns:
        obj_list : modified list with frequencies placed in sources, and 
                    passive impedances calculated
    """

    freq_set = set()
    for obj in obj_list:
        if obj.form == 'dc':
            freq_set.add(0)
    
    for line in ac_lines:
        tokens = line.split()

        if len(tokens) != 3:
            print("ERR: Invalid .ac command:", line)
            return None
        
        name = tokens[1]
        freq = get_quant(tokens[2])
        if freq == None:
            print("ERR: Frequency not a numeric type:", tokens[2])
            return None
        
        if freq not in freq_set and len(freq_set) == 0:
            freq_set.add(freq)
        elif freq not in freq_set:
            print("ERR: Multiple source frequencies not supported")
            return None
        
        for ind in range(len(obj_list)):
            obj = obj_list[ind]
            if obj.name == name:
                if obj.el_type not in 'VI' or obj.form != 'ac':
                    print("ERR: Element", name, "not an AC voltage or current source")
                    return None
                else:
                    obj.freq = freq
                    break
            
            obj_list[ind] = obj

        else:
            print("ERR: No source named", name, "defined")
            return None

        for ind in range(len(obj_list)):
            obj = obj_list[ind]
            if obj.el_type == 'L':
                obj.imp = 2*np.pi*freq*1j*obj.value
            elif obj.el_type == 'C':
                obj.imp = -1j/(2*np.pi*freq*obj.value)

            obj_list[ind] = obj
    
    return obj_list

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
        print("ERR: File not found:", path)
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
            elif not end:
                end = i
                break
        elif lines[i][:len(AC)] == AC:
            if not start:
                print("ERR: .ac found before .circuit")
                return None
            elif not end:
                print("ERR: .ac found inside .circuit and .end")
    
    
    if not start:
        print("ERR: No .circuit found")
        return None
    if not end:
        print("ERR: No .end found")
        return None

    ac_lines = []
    for i in range(end+1, len(lines)):
        if lines[i][:len(AC)] == AC:
            ac_lines.append(lines[i])
    
    
    # Get all relevant lines, reverse the order of lines
    lines = lines[start:end]

    # Parse tokens into correct fields, exit on erroneous input
    obj_list= []
    name_set = set()                        # to detect repeated names

    for line in lines:
        
        # Omit empty lines
        if line == '':
            continue

        # Get space separated tokens
        tokens = line.split()

        # Parse tokens into objects
        el = token_parse(tokens)
        
        if not el:
            return None

        elif el.name in name_set:
            print("ERR: Repeated definition for", el.name)
            return None
        
        else:
            obj_list.append(el)
    
    return [obj_list, ac_lines]

def separate_types(obj_list):
    """
    Takes list of Element objects and returns them as dictionary of 
    elements grouped by type

    Args:
        obj_list : list(Element) : List of elememts

    Returns:
        obj_dict : dict(el_type,list) : Elements grouped by type
    """

    obj_dict = {
        'R':[],
        'L':[],
        'C':[],
        'V':[],
        'I':[],
        'E':[],
        'G':[],
        'H':[],
        'F':[]
    }

    for obj in obj_list:
        obj_dict[obj.el_type].append(obj)

    return obj_dict

def get_distinct_nodes(obj_list):
    """
    Takes list of Element objects, returns dictionary of node_name:node_number pairs

    Args:
        obj_list : list(Element) : List of elements

    Returns:
        node_dict : dict(name,num) : Node name to number mapping
    """

    # Get distinct nodes using set object
    node_set = set()
    for obj in obj_list:
        node_set.add(obj.n1)
        node_set.add(obj.n2)

    if 'GND' not in node_set:
        print("ERR: No GND node found")
        return None

    # Remove GND, get node list and place GND at the beginning
    node_set.remove('GND')
    node_list = ['GND'] + sorted(list(node_set))

    # Number nodes using list index, GND gets 0
    node_dict = {node_list[i]:i for i in range(len(node_list))}

    return [node_dict, node_list]

def get_equation_matrix(obj_dict, node_dict):
    """
    Takes dictionary of objects by type and the node mapping, returns Numpy arrays representing 
    equations governing the circuit.
    WARN: Currently does not support controlled sources

    Args:
        obj_dict : dict(string, list(Element))
        node_dict : dict(string, int) : Node name to number mapping

    Returns:
        coeff_mat : N x N Numpy Array : System of equation coefficients (LHS)
        const_vec : N x 1 Numpy Array : Independent sources (RHS)
    """

    v_src_count = len(obj_dict['V'])
    v_src_map = {}
    for ind in range(len(obj_dict['V'])):
        v_src_map[obj_dict['V'][ind].name] = ind 
    
    node_count = len(node_dict.keys())-1                                   # excluded GND
    
    A = np.zeros((node_count, node_count), dtype='complex')                # top left quad
    B = np.zeros((node_count, v_src_count), dtype='complex')               # top right quad
    aux = np.zeros((v_src_count, node_count), dtype='complex')             # bottom left quad
    zero_quad = np.zeros((v_src_count,v_src_count), dtype='complex')       # bottom right quad
    const_vec_nodes = np.zeros((node_count), dtype='complex')
    const_vec_aux = np.zeros((v_src_count), dtype='complex')

    for k in obj_dict.keys():
        lst = obj_dict[k]
        for obj in lst:
            n1 = node_dict[obj.n1]-1
            n2 = node_dict[obj.n2]-1
            value = obj.value

            if obj.el_type in 'RLC':
                if obj.imp == None:   # dc
                    if obj.el_type == 'L':
                        value = EPS
                    elif obj.el_type == 'C':
                        value = INF
                else:
                    value = obj.imp + EPS

                if n1 != -1 and n2 != -1:
                    A[n1,n1] += 1/value
                    A[n1,n2] += -1/value
                    A[n2,n1] += -1/value
                    A[n2,n2] += 1/value

                elif n1 == -1:
                    A[n2,n2] += 1/value

                elif n2 == -1:
                    A[n1,n1] += 1/value

            elif obj.el_type == 'V':
                nv = v_src_map[obj.name]
                if n1 != -1 and n2 != -1:
                    B[n1,nv] += 1
                    B[n2,nv] += -1
                    aux[nv,n1] += 1
                    aux[nv,n2] += -1
                
                elif n1 == -1:
                    B[n2,nv] += -1
                    aux[nv,n2] += -1
                
                elif n2 == -1:
                    B[n1,nv] += 1
                    aux[nv,n1] += 1

                const_vec_aux[nv] += value

            elif obj.el_type == 'I':
                if n1 != -1 and n2 != -1:
                    const_vec_nodes[n1] += obj.value
                    const_vec_nodes[n2] += -obj.value
                elif n1 == -1:
                    const_vec_nodes[n2] += -obj.value
                elif n2 == -1:
                    const_vec_nodes[n1] += obj.value

            elif obj.el_type in 'EGHF':
                print("ERR: Controlled sources not supported yet")
                return None

    top_half = np.concatenate((A, B), axis=1)  
    bot_half = np.concatenate((aux, zero_quad), axis=1)
    coeff_mat = np.concatenate((top_half, bot_half), axis=0)
    const_vec = np.concatenate((const_vec_nodes, const_vec_aux), axis=0)

    return [coeff_mat, const_vec]

def solve_matrix(M, b):
    """
    Takes M and b in a linear system Mx=b, and returns x = M_inv * b

    Args:
        M : Q x Q Numpy Array : Coefficient matrix
        b : Q x 1 Numpy Array : Constant vector

    Returns:
        x : Q x 1 Numpy Array : Variable vector
    """

    try:
        x = np.linalg.solve(M, b)
    except np.LinAlgError:
        print("ERR: Matrix is singular")
        return None

    if not np.allclose(np.dot(M, x), b):
        print("ERR: Matrix is inconsistent (most likely with the independent sources)")
        return None
    
    return x

def display_sol(sol, node_list, v_src_list):
    """
    Takes solution vector and variable names in the form of node voltages and voltage
    source currents, and prints them in user-friendly fashion

    Args:
        sol : 1D Numpy Array : Solution vector
        node_list : list(string) : List of nodes
        v_src_lst : list(Element) : List of voltage sources

    Returns:
        None
    """

    print("\n")
    node_list = node_list[1:]
    for i in range(len(node_list)):
        print("V_" + node_list[i] + ": ", ffs(np.real(sol[i]), precision=5), '+', ffs(np.imag(sol[i]), precision=5)+'j')

    for i in range(len(v_src_list)):
        v = v_src_list[i]
        print("I_" + v.name + ": ", ffs(np.real(sol[len(node_list)+i]), precision=5), '+', ffs(np.imag(sol[len(node_list)+i]), precision=5)+'j')
    print("\n")

def EXIT_ON_NONE(var):
    """
    If var is None, it means an error has happened (and been reported by earlier prints), 
    so exit execution with code 1

    Args:
        var : Variable to check

    Returns:
        None
    """

    try:
        if not var.any():
            pass
    except:
        if var == None:
            sys.exit(1)

if __name__ == '__main__':
    arg_count = len(sys.argv)
    if arg_count > 2:
        print("ERR: Too many arguments (expected 1 input file)")
        sys.exit(1)
    elif arg_count < 2:
        print("ERR: Too less arguments (expected 1 input file)")
        sys.exit(1)
    
    ret = netlist_to_data(sys.argv[1])
    EXIT_ON_NONE(ret)
    obj_list, ac_lines = ret
    EXIT_ON_NONE(obj_list)
    EXIT_ON_NONE(ac_lines)
    print("DEBUG: File reading done")
    
    obj_list = ac_parse(obj_list, ac_lines)
    EXIT_ON_NONE(obj_list)
    print("DEBUG: AC lines parsed")
    
    obj_dict = separate_types(obj_list)
    EXIT_ON_NONE(obj_dict)
    print("DEBUG: Elements segregated")
    
    node_dict, node_list = get_distinct_nodes(obj_list)
    EXIT_ON_NONE(node_dict)
    print("DEBUG: Distinct nodes obtained")
    
    ret = get_equation_matrix(obj_dict, node_dict)
    EXIT_ON_NONE(ret)
    coeff_mat, const_vec = ret
    EXIT_ON_NONE(coeff_mat)
    EXIT_ON_NONE(const_vec)
    print("DEBUG: Equation matrix obtained")

    sol = solve_matrix(coeff_mat, const_vec)
    EXIT_ON_NONE(sol)
    print("DEBUG: Solution obtained")

    display_sol(sol, node_list, obj_dict['V'])
    print('DEBUG: Done')

    # normal exit code 0
    sys.exit(0)
    