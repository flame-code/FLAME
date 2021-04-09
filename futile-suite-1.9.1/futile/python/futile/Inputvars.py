"""Handle the input variable specifications.

This module is the python complement to the specification of the input variables of a 
file as provided by the :f:mod:`f_input_file` module. It uses the same syntax as defined there 
and make possible the interplay between a python-based serialization of the input dictionary 
and the construction of 

"""

class InputVariable():
    """
    Define a input variable of a library.
    Such object can be initialized and inspected from the dictionary used by futile 
    and then initialized according to the provided specification.
    """
    def __init__(self,name,spec):
        """
        Define a specification of the input variable, according to the conventions
        employed in the :f:mod:`f_input_file` module.
        
        Args:
           name (str): the name of the input variable
           spec (dict): dictionary describing the specification of the variables

        Todo: 
           Extensively specify the conventions somewhere

        """
        from futile.Utils import kw_pop
        if 'default' not in spec or 'COMMENT' not in spec:
            raise ValueError("Invalid specification of the variable '"+name+"'")
        self.name=name
        self.spec,self.valid_range=kw_pop('RANGE',None,**spec)
        self.spec,self.valid_values=kw_pop('EXCLUSIVE',None,**self.spec)
        self.spec,desc=kw_pop('DESCRIPTION','None',**self.spec)
        self.spec,self.profile_from=kw_pop('PROFILE_FROM',None,**self.spec)
        self.spec,self.condition=kw_pop('CONDITION',None,**self.spec)
        self.impose_profile_to=[]
        self._val_from_profile('default')
        self.__doc__=self.name+":\n"+self.spec.pop("COMMENT")+"\nDescription:\n"+desc
    def _val(self,val):
        if self.valid_range and             (val < self.valid_range[0] or val > self.valid_range[1]):
                raise ValueError("Variable '"+ str(val)+ "' not in allowed range " + str(self.valid_range))
        if self.valid_values and (val not in self.valid_values):
            raise ValueError("Variable '" + str(val) + "' not in allowed values " + str(self.valid_values))
        self.value=val
    def _val_from_profile(self,profile):
        self.profile=profile
        val=self.spec[profile]
        if type(val) == str and val in self.spec: val=self.spec[val]
        self._val(val)
        #if the profile chosen is also present in the input variables that 
        #depends from that then impose it
        for name,var in self.impose_profile_to:
            if self.profile in var.spec: var.set(self.profile)
    def is_valid(self):
        if self.condition is not None and "VARIABLE" in self.condition:   
            #raise ValueError("Variable '"+self.name+"' requires to be bound with master_key '"+\
            #                 self.condition["MASTER_KEY"]+"';\n use the set_master_variable method")
            var=self.condition["VARIABLE"]
            when=self.condition.get("WHEN")
            if when and (var.profile not in when and var.value not in when):
                return False #raise ValueError("Variable '"+self.name+"' cannot be set as master key '"+\ self.condition["MASTER_KEY"]+"' does not have the allowed values")
            when_not=self.condition.get("WHEN_NOT")
            if when_not and (var.profile in when_not or var.value in when_not):
                return False #raise ValueError("Variable '"+self.name+"' cannot be set as master key '"+self.condition["MASTER_KEY"]+"' has incompatible values")
        return True
    def _val_from_user(self,val):
        self.profile='__USER__'
        #search if some profile corresponds to a value
        for profile,value in self.spec.items():
            if value == val:
                self._val_from_profile(profile)
                return
        self._val(val)
    def set(self,val):
        """Set the value of the input variable"""
        if val in self.spec:
            self._val_from_profile(profile=val)
        else:
            self._val_from_user(val)
    def set_dependent_variable(self,var):
        """Set the dependent variable from which impose the profile"""
        self.impose_profile_to.append((var.name,var))
    def set_master_variable(self,var):
        """Set the variable which is activated when the present has suitable values or profiles"""
        if self.condition: 
            if var.name != self.condition["MASTER_KEY"]:
                raise ValueError("Variable '"+self.name+"' requires to be bound with '"+                                self.condition["MASTER_KEY"]+"', whereas '"+                                var.name+"' has been provided")
            self.condition["VARIABLE"]=var
    def __repr__(self):
        if self.profile != '__USER__':
            return self.profile
        else:
            return str(self.value)

