function out=binary(idel)

ON=max(idel);
OFF=min(idel);

length_t= size(idel,1);

bin=zeros(length_t,1);

for i=1:length_t
    if idel(i,1)==ON
        bin(i,1)=1;
        if idel(i,1)==OFF
           bin(i,1)=0; 
        end
    end
end


out=bin;
end


