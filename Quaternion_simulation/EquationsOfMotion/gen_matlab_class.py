import sympy as sp

def MatlabFunction(function,fun_name,assignto,coordinates,speeds,inputs,parameters):
    
    text_file = open(f"{fun_name}.m","w")
    with open(f"{fun_name}.m","w") as text_file:
        # function .. = .. ()
        header = f"function {assignto} = {fun_name}(q,u,inputs,model)" #,u,act,model,opt_var)
        
        print(header,file=text_file)
        
        for par_name, par_tuple in parameters.items():
            if not type(par_tuple) == list:
                print(f'{par_name} = model.{par_name};',file=text_file)
            else:
                for i, elem in enumerate(var_tuple):
                    str_res = f'{elem} = model.{par_name}(:,{i+1});'
                    print(str_res,file=text_file)
        
        for i, coord in enumerate(coordinates):
            print(f"{str(coord)} = q({i+1},:);", file=text_file)
        for i, speed in enumerate(speeds):
            print(f"{str(speed)} = u({i+1},:);", file=text_file)
                    
        for i, act in enumerate(inputs):
            print(f"{str(act)} = inputs({i+1},:);", file=text_file)

                
        sub_exprs, simplified_rhs = sp.cse(function)
        for var, expr in sub_exprs:
            print('%s = %s;' % (sp.octave_code(var),sp.octave_code(expr)),file=text_file)
        print('%s' % sp.octave_code(sp.Matrix([simplified_rhs]), assign_to = assignto),file=text_file)