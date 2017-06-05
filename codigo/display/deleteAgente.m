function [ListBox,Array] = deleteAgente(LB,A)
	ListBox = LB;
	Array = A;
	if ListBox.Value == 1 && length(ListBox.Items) == 1
		ListBox.Items = {};
		ListBox.ItemsData = [];
		Array = [];
		ListBox.Value = {};
	else
		ind = setdiff(1:length(ListBox.Items),ListBox.Value);
		ListBox.Items = ListBox.Items(ind);
		ListBox.ItemsData = (1:length(ListBox.Items));
		Array = Array(ind);
		ListBox.Value = 1;
	end
end
