function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'aod', 'james_2pi');
end
obj = schemaObject;
