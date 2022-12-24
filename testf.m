 File=load('matlab');
 D_Kp={};
 D_Kp_NN =[];
for i=1:1471

 if File.Data(i,3)== 1
  D_Kp(i)={'ZE'};
  D_Kp_NN(i,:)=[1,0,0,0,0];
 elseif File.Data(i,3)== 2
   D_Kp(i)={'S'};
  D_Kp_NN(i,:)=[0,1,0,0,0];
  elseif File.Data(i,3)== 3
   D_Kp(i)={'M'};
 D_Kp_NN(i,:)=[0,0,1,0,0];
  elseif File.Data(i,3)== 4
   D_Kp(i)={'B'};
  D_Kp_NN(i,:)=[0,0,0,1,0];
   elseif File.Data(i,3)== 5
   D_Kp(i)={'VB'};
  D_Kp_NN(i,:)=[0,0,0,0,1];
 end
end
  D_Ki={};
 D_Ki_NN =[];
for i=1:1471

 if File.Data(i,4)== 1
   D_Ki(i)={'ZE'};
  D_Ki_NN(i,:)=[1,0,0,0,0];
 elseif File.Data(i,4)== 2
    D_Ki(i)={'S'};
  D_Ki_NN(i,:)=[0,1,0,0,0];
  elseif File.Data(i,4)== 3
    D_Ki(i)={'M'};
 D_Ki_NN(i,:)=[0,0,1,0,0];
  elseif  File.Data(i,4)== 4
    D_Ki(i)={'B'};
  D_Ki_NN(i,:)=[0,0,0,1,0];
   elseif File.Data(i,4)== 5
    D_Ki(i)={'VB'};
  D_Ki_NN(i,:)=[0,0,0,0,1];
 end
end
  D_Kd={};
 D_Kd_NN =[];
for i=1:1471

 if File.Data(i,5)== 1
   D_Kd(i)={'ZE'};
  D_Kd_NN(i,:)=[1,0,0,0,0];
 elseif  File.Data(i,5)== 2
    D_Kd(i)={'S'};
  D_Kd_NN(i,:)=[0,1,0,0,0];
  elseif File.Data(i,5)== 3
    D_Kd(i)={'M'};
 D_Kd_NN(i,:)=[0,0,1,0,0];
  elseif File.Data(i,5)== 4
    D_Kd(i)={'B'};
  D_Kd_NN(i,:)=[0,0,0,1,0];
   elseif File.Data(i,5)== 5
    D_Kd(i)={'VB'};
  D_Kd_NN(i,:)=[0,0,0,0,1];
 end
end
D_Labels= [D_Kp; D_Ki; D_Kd;]
D_NN_Labels=[D_Kp_NN_Labels; D_Ki_NN; D_Kd_NN];