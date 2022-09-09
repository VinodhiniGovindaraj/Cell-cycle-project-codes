% Paper: Transcriptional fluctuations govern the serum dependent cell cycle duration heterogeneities in Mammalian cells
% Author: Vinodhini Govindaraj, Subrot Sarma, Atharva Karulkar, Rahul Purwar and Sandip Kar
% e-mail about the code: vinodhinigovindaraj@gmail.com,sandipkar@iitb.ac.in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [akn] = trvariation(akm,ku,vf)
%TRVARIATION Summary of this function goes here
%   Detailed explanation goes here

	akn=akm;
    
    
for j=1:ku

 		rk3 = rand();

        if (rk3 <= 0.5)
            if (rk3>=0.4)
			akn(j)=akn(j)-(0.01*vf*akn(j));
			
            elseif (rk3>=0.3)
			akn(j)=akn(j)-(0.02*vf*akn(j));
			

            elseif (rk3>=0.2)
			akn(j)=akn(j)-(0.03*vf*akn(j));

            elseif (rk3>=0.1)
			akn(j)=akn(j)-(0.04*vf*akn(j));
			
            elseif (rk3>=0.0)
			akn(j)=akn(j)-(0.05*vf*akn(j));
			
            end
            
        elseif (rk3 >= 0.5)

            if (rk3<=0.6)
			akn(j)=akn(j)+(0.01*vf*akn(j));
			
            elseif (rk3<=0.7)
			akn(j)=akn(j)+(0.02*vf*akn(j));
		

            elseif (rk3<=0.8)
			akn(j)=akn(j)+(0.03*vf*akn(j));
		
            elseif (rk3<=0.9)
			akn(j)=akn(j)+(0.04*vf*akn(j));
		
            elseif (rk3<=1.0)
			akn(j)=akn(j)+(0.05*vf*akn(j));
		
            end 

        end 
end



