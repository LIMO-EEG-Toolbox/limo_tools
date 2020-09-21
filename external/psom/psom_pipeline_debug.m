function psom_pipeline_debug(commands)

fields1 = fieldnames(commands);
for iField1 = 1:length(fields1)
    fprintf('Executing LIMO batch pipeline %s ...\n', fields1{iField1});
    
    fields2 = setdiff( fieldnames(commands.(fields1{iField1})), 'command');
    for iField2 = 1:length(fields2)
        eval( [ fields2{iField2} '= commands.(fields1{iField1}).(fields2{iField2});' ] );
    end
    eval( [ commands.(fields1{iField1}).command ] );
end