function [horas] = loadHoras(Config)

horas = cell(Config.Etapas,1);
hh_num = Config.Etapas/4;
hh = reshape(strrep([' ' num2str(mod((0:Config.Etapas/4-1),24),'%2u')],' ','0'),2,hh_num)';
mm_num = 4;
mm = reshape(strrep([' ' num2str((0:15:45),'%2d')], ' ', '0'),2,4)';

for i = 1:Config.Etapas
	hh_i = idivide(int32(i-1),4)+1;
	mm_i = mod(i,mm_num);
	if mm_i == 0;
		mm_i = 4;
	end
	horas{i} = [hh(hh_i,:) ':' mm(mm_i,:)];
end