function h = myplot (obj)
  switch obj.Type
    case "oneway"
      h = boxplot (obj.Data, obj.Groups);

    case "twoway"
      % Later: add interaction plots for anova2

    case "nway"
      % Later: add factorial plots for anovan
  end
end
