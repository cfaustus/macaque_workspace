

%throughout this model, each alpha globin cluster is represented as two
%digits, representing HBA1 and HBA2.

%The possible alpha globin types at HBA1 and HBA2 are types a, b and c (see
%main text). Throughout this model, type a is indicated by the number 1,
%type b by the number 2 and type c by the number 3.

%A complete macaque genotype is given by 4 numbers -the first two
%representing HBA1 and HBA2 on one of the macaque's chromosomes ,the second
%two representing HBA1 and HBA2 on the other chromosome.

%So, e.g. a macaque encoding type a globin in all its alpha globin genes is
%represented by "1111".

%In figure 3A-E, the population started with all wild type alpha globin:
initpop=ones(10000,4);

%The following code (currently commented out) provides the starting population
%with a pattern of HBA1 and HBA2 mutations (as used in figure 3F and 3G).

% for i=1:1:10000
%     
%  if (unifrnd(0,1)<0.15)
%      
%      check=unifrnd(0,1);
%      
%      if check<(1/3)
%      
%      initpop(i,1)=1;
%      initpop(i,2)=2;
%      
%      end
%      
%      if check>=(1/3) && check<(2/3)
%      
%      initpop(i,1)=3;
%      initpop(i,2)=1;
%      
%      end
%      
%      
%      if check>=(2/3) 
%          
%      initpop(i,1)=3;
%      initpop(i,2)=2;
%      
%      end
%      
%      
%  end
%  
%  
%  
%   if (unifrnd(0,1)<0.15)
%      
%      check=unifrnd(0,1);
%      
%      if check<(1/3)
%      
%      initpop(i,3)=1;
%      initpop(i,4)=2;
%      
%      end
%      
%      if check>=(1/3) && check<(2/3)
%      
%      initpop(i,3)=3;
%      initpop(i,4)=1;
%      
%      end
%      
%      
%      if check>=(2/3) 
%          
%      initpop(i,3)=3;
%      initpop(i,4)=2;
%      
%      end
%      
%      
%  end
%      
% end


%for this part of the code to work, "AlphaGlobinPopulationGeneticModel.c"
%must have been compiled as a MEX file. AlphaGlobinPopulationGeneticModel
%takes as one of its input arguments the starting population conditions,
%but as a vector, so "startpop" needs to be reshaped. Likewise, the outputs
%of AlphaGlobinPopulationGeneticModel are vectors which need to be reshaped
%(lines ??) to make them easier to interpret.

startpop=reshape(initpop,1,40000);

%the other imput for AlphaGlobinPopulationGeneticModel are the following
%parameters, which correspond to the parameters in the Methods as shown here: 

% RunTime= t
% MutR = de novo mutation rate (never used, always = 0 in this paper).
% MalSel= d
% prot = p
% costAbnormal = k
% recombR= r
% convR=c
% migR= m



%another input  for AlphaGlobinPopulationGeneticModel is a random
%number, which will later be added to the current time in order to provide 
%a random seed for the random number generator in the population genetic model. 
randomnumber=unifrnd(1,100000);

parameters= [RunTime,MutR,MalSel,prot,costAbnormal,recombR,convR, migR,randomnumber];


[resline,haplos]=AlphaGlobinPopulationGeneticModel(parameters,startpop);

ResultPop= reshape(resline,4,10000)';
haplotimeline= reshape(haplos,9,RunTime/100)';

%"ResultPop" contains 10000 rows of four numbers. Each row represents an
%individual macaque genotype (see notes at beginning of this script).

%The 9 columns of "haplotimeline" contains the frequencies of the 9
%possible alpha globin clusters, measured every 100 generations of the simulation.

%The order of the 9 alpha globin clusters are:

% 11  (a a)
% 12  (a b)
% 13  (a c)
% 21  (b a)
% 22  (b b)
% 23  (b c)
% 31  (c a)
% 32  (c b)
% 33  (c c)



