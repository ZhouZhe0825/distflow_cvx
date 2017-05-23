function copyrename(file,renfile,newdir)
    copyfile(file, newdir);
    spl = strsplit(file, '\');
    if ~exist([newdir,'\',renfile],'file')
        movefile([newdir,'\',spl{end}],[newdir,'\',renfile]);
    end
end
