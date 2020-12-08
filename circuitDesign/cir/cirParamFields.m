function cellArray = cirParamFields(circuit, eltype)
switch eltype
  case {'nmos', 'pmos'}
    cellArray = circuit.defaults.(eltype).fieldValues{4};
  case {'res', 'cap', 'ind', 'v', 'i'}
    cellArray = circuit.defaults.(eltype).fieldValues{2};
otherwise
    error('no defaults exist for element type %s', eltype);
end

    
