function measOut = pullMeasFromList(obsList,type)

measOut = [];

for idx = 1:length(obsList)
    obsi = obsList{idx};
    if obsi.type == type
        measOut = obsi;
    end
end


end