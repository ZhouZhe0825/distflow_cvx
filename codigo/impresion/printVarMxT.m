function printVarMxT(A, nodos, Header, outFilename, sheetName)
	
	sheetData = [Header;[nodos num2cell(A)]];
	xlswrite([outFilename, '.xlsx'], sheetData, sheetName);

end
