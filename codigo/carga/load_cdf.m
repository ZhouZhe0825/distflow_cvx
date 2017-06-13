function [Bus, Branch] = load_cdf(filename);

%% Sample program to read data from IEEE Common Data Format (Tried on 30 Bus Data)
%% Bus Data and Line Data are read into matrices which can be used for load flow
%% Please Refer: http://www.ee.washington.edu/research/pstca/pf30/pg_tca30bus.htm
%% Coded by: Krishnanand K.R., ECE, National University of Singapore
%% Supervisor: Prof. Sanjib K. Panda, ECE, National University of Singapore
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fid = fopen('power_system_test\ieee14cdf.txt');
% fid = fopen('power_system_test\ieee30cdf.txt');
% fid = fopen('power_system_test\ieee57cdf.txt');
% fid = fopen('power_system_test\ieee118cdf.txt');

fid = fopen(filename);

%% BUS Data
Line_String_Complete=fgetl(fid);
Line_String_Complete=fgetl(fid);
Bus.Num = [];
Bus.Name = [];
Bus.Lfan = [];
Bus.Lzn = [];
Bus.Type = [];
Bus.Fv = [];
Bus.Fa = [];
Bus.Lmw = [];
Bus.Lmvar = [];
Bus.Gmw = [];
Bus.Gmvar = [];
Bus.Bkv = [];
Bus.Dv = [];
Bus.Maxmvar = [];
Bus.Minmvar = [];
Bus.Scg = [];
Bus.Scb = [];
Bus.Rcbn = [];
No_of_Buses=0;
while ischar(Line_String_Complete)
	Line_String_Complete=fgetl(fid);
	%     disp(['#' Line_String_Complete '#'] );
	if(strcmp(Line_String_Complete(1:4),'-999')==1)
		break;
	end
	Bus.Num = [Bus.Num; str2num(Line_String_Complete(1:4))];
	Bus.Name = [Bus.Name; Line_String_Complete(6:17)];
	Bus.Lfan = [Bus.Lfan; str2num(Line_String_Complete(19:20))];
	Bus.Lzn = [Bus.Lzn; str2num(Line_String_Complete(21:23))];
	Bus.Type = [Bus.Type; str2num(Line_String_Complete(25:26))];
	Bus.Fv = [Bus.Fv; str2num(Line_String_Complete(28:33))];
	Bus.Fa = [Bus.Fa; str2num(Line_String_Complete(34:40))];
	Bus.Lmw = [Bus.Lmw; str2num(Line_String_Complete(41:49))];
	Bus.Lmvar = [Bus.Lmvar; str2num(Line_String_Complete(50:59))];
	Bus.Gmw = [Bus.Gmw; str2num(Line_String_Complete(60:67))];
	Bus.Gmvar = [Bus.Gmvar; str2num(Line_String_Complete(68:75))];
	Bus.Bkv = [Bus.Bkv; str2num(Line_String_Complete(77:83))];
	Bus.Dv = [Bus.Dv; str2num(Line_String_Complete(85:90))];
	Bus.Maxmvar = [Bus.Maxmvar; str2num(Line_String_Complete(91:98))];
	Bus.Minmvar = [Bus.Minmvar; str2num(Line_String_Complete(99:106))];
	Bus.Scg = [Bus.Scg; str2num(Line_String_Complete(107:114))];
	Bus.Scb = [Bus.Scb; str2num(Line_String_Complete(115:122))];
	Bus.Rcbn = [Bus.Rcbn; str2num(Line_String_Complete(124:end))];

	No_of_Buses=No_of_Buses+1;
end

%% Line Data
Line_String_Complete=fgetl(fid);
Branch.Tbn = [];
Branch.Zbn = [];
Branch.Lfa = [];
Branch.lz = [];
Branch.Circuit = [];
Branch.Type = [];
Branch.Brr = [];
Branch.Brx = [];
Branch.Lcb = [];
Branch.Lmva1 = [];
Branch.Lmva2 = [];
Branch.Lmva3 = [];
Branch.Cbn = [];
Branch.Side = [];
Branch.Tftr = [];
Branch.Tfa = [];
Branch.Mint = [];
Branch.Maxt = [];
Branch.Ss = [];
Branch.Minv = [];
Branch.Maxv = [];
No_of_Branchs=0;
while ischar(Line_String_Complete)
	Line_String_Complete=fgetl(fid);
	%     disp(['#' Line_String_Complete '#'] );
	if(strcmp(Line_String_Complete(1:4),'-999')==1)
		break;
	end
	Branch.Tbn = [Branch.Tbn; str2num(Line_String_Complete(1:4))];
	Branch.Zbn = [Branch.Zbn; str2num(Line_String_Complete(6:9))];
	Branch.Lfa = [Branch.Lfa; str2num(Line_String_Complete(11:12))];
	Branch.lz = [Branch.lz; str2num(Line_String_Complete(13:15))];
	Branch.Circuit = [Branch.Circuit; str2num(Line_String_Complete(17))];
	Branch.Type = [Branch.Type; str2num(Line_String_Complete(19))];
	Branch.Brr = [Branch.Brr; str2num(Line_String_Complete(20:29))];
	Branch.Brx = [Branch.Brx; str2num(Line_String_Complete(30:40))];
	Branch.Lcb = [Branch.Lcb; str2num(Line_String_Complete(41:50))];
	Branch.Lmva1 = [Branch.Lmva1; str2num(Line_String_Complete(51:55))];
	Branch.Lmva2 = [Branch.Lmva2; str2num(Line_String_Complete(57:61))];
	Branch.Lmva3 = [Branch.Lmva3; str2num(Line_String_Complete(63:67))];
	Branch.Cbn = [Branch.Cbn; str2num(Line_String_Complete(69:72))];
	Branch.Side = [Branch.Side; str2num(Line_String_Complete(74))];
	Branch.Tftr = [Branch.Tftr; str2num(Line_String_Complete(77:82))];
	Branch.Tfa = [Branch.Tfa; str2num(Line_String_Complete(84:90))];
	Branch.Mint = [Branch.Mint; str2num(Line_String_Complete(91:97))];
	Branch.Maxt = [Branch.Maxt; str2num(Line_String_Complete(98:104))];
	Branch.Ss = [Branch.Ss; str2num(Line_String_Complete(106:111))];
	Branch.Minv = [Branch.Minv; str2num(Line_String_Complete(113:119))];
	Branch.Maxv = [Branch.Maxv; str2num(Line_String_Complete(120:end))];
    No_of_Branchs=No_of_Branchs+1;
end

fclose(fid);


