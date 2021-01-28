function mosTable = tablePreprocess(mosTable, wmin, lmin, wref, wcrit, ...
    vddMax, modelname, modelStruct, techName)
%TABLEPREPROCESS preprocessing of a table of parameters of a MOS
%transistor.
%
%  TABLEPREPROCESS(MOSTABLE, WMIN, LMIN, WREF, WCRIT, VDDMAX, MODELNAME,
%  MODELSTRUCT, TECHNAME) changes a given raw table MOSTABLE, such that it is
%  uniform for every technology. For example, the generated raw table can be 
%  different depending on the MOS model: for a BSIM model, Spectre
%  simulations generate different operating point parameters 
%  (and sometimes with different names), than for a MM11 model.
%  
%  The following fields are added to the MOSTABLE
%    - Input
%    - Model. This contains several model parameters, namely the ones that 
%      are needed for computations in existing functions, such as tox, VTO
%      (if existing, does not exist in MM11), all parameters necessary
%      to calculate junction capacitors for bulk CMOS, ...
%      The model parameters that will be put in this field are the ones
%      that correspond to a field of the structure MODELSTRUCT that is
%      given as an argument to this function. For example, if MODELSTRUCT
%      has a field MODELSTRUCT.tox with value 2e-9, then a field 
%      MOSTABLE.Model.tox will be created and it has a value of 2e-9.
%    - extra fields of MOSTABLE.Info: minimum allowed channel length in meters 
%      and width in meters (= arguments LMIN and WMIN), reference width in
%      meters (argument WREF), critical width in meters (argument WCRIT), 
%      maximum allowed power supply voltage in Volts (argument VDDMAX),
%      name of the model (argument MODELNAME, e.g. 'MOS Model 11') and name
%      of the technology (argument TECHNAME, e.g. 'IMEC CMOS 90 nm')
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

mosTable.Info.tech.name = techName;
mosTable.Info.inputs = {'lg', 'vgs', 'vds', 'vsb'};


mosTable.Input.vsb.step=-1*mosTable.Info.VbsStep;
mosTable.Input.vsb.final=-1*mosTable.Info.VbsFinal;
mosTable.Input.vsb.init = 0;
mosTable.Input.vsb.uniformStep = 1;
mosTable.Input.vsb.array = 0:mosTable.Input.vsb.step:mosTable.Input.vsb.final;
mosTable.Input.vsb.nInputs = length(mosTable.Input.vsb.array);

mosTable.Input.vgs.step = mosTable.Info.VgsStep;
mosTable.Input.vgs.final = mosTable.Info.VgsFinal;
mosTable.Input.vgs.init = 0;
mosTable.Input.vgs.uniformStep = 1;
mosTable.Input.vgs.array = 0:mosTable.Info.VgsStep:mosTable.Info.VgsFinal;
mosTable.Input.vgs.nInputs = length(mosTable.Input.vgs.array);

mosTable.Input.vds.step = mosTable.Info.VdsStep;
mosTable.Input.vds.final = mosTable.Info.VdsFinal;
mosTable.Input.vds.init = 0;
mosTable.Input.vds.uniformStep = 1;
mosTable.Input.vds.array = 0:mosTable.Info.VdsStep:mosTable.Info.VdsFinal;
mosTable.Input.vds.nInputs = length(mosTable.Input.vds.array);

mosTable.Input.lg.step = NaN;
mosTable.Input.lg.final = 1e-6*mosTable.Info.Larray(length(mosTable.Info.Larray));
mosTable.Input.lg.init = 1e-6*mosTable.Info.Larray(1);
mosTable.Input.lg.uniformStep = 0;
mosTable.Input.lg.array = 1e-6*mosTable.Info.Larray;
mosTable.Input.lg.nInputs = length(mosTable.Input.lg.array);

% modification of mosTable.Table:
mosTable.Table.cdbI = mosTable.Table.cdb;
mosTable.Table.csbI = mosTable.Table.csb;

mosTable.Table = rmfield(mosTable.Table, {'vgs', 'vds', 'vbs', 'cdb', 'csb'});

% filling the field mosTable.Model:
names = fieldnames(modelStruct);
for i = 1:length(names)
  mosTable.Model.(char(names(i))) = modelStruct.(char(names(i)));
end

% modification of mosTable.Info:
mosTable.Info.spectreModelFile = mosTable.Info.GenericModel;
mosTable.Info.vddMax = vddMax;
mosTable.Info.WIndepParams = {'vth', 'vdsat'};
mosTable.Info = rmfield(mosTable.Info, {'VbsStep', 'VbsFinal', 'VgsStep', 'VgsFinal', ...
    'VdsStep', 'VdsFinal', 'Larray', 'GenericModel'});
mosTable.Info.modelName = modelname;
mosTable.Info.modelParams = fieldnames(mosTable.Model);
mosTable.Info.Wmin = wmin;
mosTable.Info.Lmin = lmin;
mosTable.Info.Wref = wref;
mosTable.Info.Wcrit = wcrit;
