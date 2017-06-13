function [pPvg] = calculoPotenciaSolar(temp, pPvg_low, pPvg_top)

	pPvg = (temp - min(temp))/max(temp)*pPvg_top + pPvg_low; % potencia de generacion del solar, por sale desde la temperatura

end