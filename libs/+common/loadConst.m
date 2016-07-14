function ret = loadConst(group)
  % litrconst: return literal constants
  % litr = litrconst().ja;
  % litr.ordinal{1}
  if nargin == 0, group = ''; end
  switch group
    case 'env'
      ret = ReadYaml('res/envconst.yaml');
    case 'exp'
      ret = ReadYaml('res/expconst.yaml');
    case 'litr'
      ret = ReadYaml('res/litrconst.yaml');
    case 'data'
      ret = ReadYaml('res/dataconst.yaml');
    otherwise
      ret = ReadYaml('res/litrconst.yaml');
  end
end
