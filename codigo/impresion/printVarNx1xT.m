function printVarNx1xT(A, nodos, Header, outFilename, sheetName)
	
	sheetData = [Header;[nodos num2cell(squeeze(A))]];
	xlswrite([outFilename, '.xlsx'], sheetData, sheetName);

end
