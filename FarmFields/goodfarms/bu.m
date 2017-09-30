clear
clc
a=dir;
z=1;
for i=1:length(a)
    n=a(i).name;
    if(length(n)>5)
        load(n)
        if(z<10)
        str=['good000' num2str(z)];
        else
        str=['good00' num2str(z)];
        end
        save(str)
        z=z+1;
    end
end

    