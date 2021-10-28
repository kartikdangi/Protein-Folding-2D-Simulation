%In this code, protein folding simulation is done with grid crossing and displayed dynamically after each 2000 monte carlo moves by 6 polymers
%Dynamic display is shown after each 2000 moves so as to reduce plotting load on matlab due to large number of iterations and to reduce time taken by matlab.
%Note : Use Energy_A,Energy_B,Energy_C functions in Line 139 with different Interaction energy value in Line 12 to compute structures for Part A,B,C accordingly.
close all;clc;%Code takes around 7/8 minutes for 10^6 iterations with dynamic display after each 2000 iterations.
SqDim = 40;
%Defines the dimension of the square lattice

%Notes:ORANGE ERROR shown is due to unused functions Energy_B and Energy_C which can be used as per user need.

StableEnergyIndex = [];
StableEnergyIndex(1) = 0;
%StableEnergyIndex stores all iterations when it satisfies metropolis criteria and lattice conformation changes.

InteractionEnergy = -3;
%It stores the value of interaction energy between two residues.

N = 1000000;%Can also be written 10^6
%N stores the number of Iterations 

EInitial = 0;
%It stores the Inital energy of lattice system when no residue interacts, which stands as the relative basis for further calculations of energy for the lattice.

Lattice(:,:,1) = [ 3 4;3 5;3 6;3 7;];
Lattice(:,:,2) = [ 15 23;16 23;17 23;18 23;];
Lattice(:,:,3) = [ 20 15;21 15;22 15;23 15;];
Lattice(:,:,4) = [ 15 17;16 17;17 17;18 17;];
Lattice(:,:,5) = [ 10 5;11 5;12 5;13 5;];
Lattice(:,:,6) = [ 5 12;5 13;5 14;5 15;];
%Above assignment creates a 3D array, with each page (Total - 6 pages for 6 proteins) of the array containing the locations of four residues for that specific protein in form of 4Rows X 2Cols, where each Row denotes one residue of that protein with 1st and 2nd column storing the X and Y coordinates for the residue on 2D grid 

InitialSystemEnergy = [];InitialSystemEnergy(1) = EInitial;
%Above array stores the energy of the lattice system and continues to append the Energy of lattice whenever Protein Movements take place.
figure(1);set(gcf, 'Position',  [50, 50, 800, 700]); xh = get(gca,'xlabel') ; p = get(xh,'position');p(2) = p(2)-1 ; set(xh,'position',p);
for Iteration = 1:N
    Lattice2 = Lattice(:,:,:);
    r = randi(2,1,1);%This random number is useful in looping over protein, without giving any priority to any polymer. 
    if r == 1
       m = 1;n=1;o=6;
    else 
      m = 6;n=-1;o=1;
    end    
    %A copy of original lattice is made after each iteration in order to make all moves
    %considering it as the orignial
    for protein_num = m:n:o
    %Above for loops over each protein for movement in each iteration from 1 to N.    
        Translate = randi(2,1,1);
        %If Translate == 1, we won't do translation motion for our selected protein
        %If Translate == 2, we will do translation motion for our selected protein
        TypeOfMove = randi(2,1,1);
        if Translate == 1
            %If Translate = 2,selected protein will try to do a translatinal move.
            p = -3+2*randi(2,1,1);
            %Above formula gives either p = -1 or p = 1            
            %p defines the direction our protein may try to move randomly either +x/-x/+y/-y in case no other protein is already present there.               
            checkOccupancy = 0;
            directionXorY = randi(2,1,1);
            %directionXorY selects randomly x or y direction
            if directionXorY == 1%Selected protein will translate along x axis/horizontally
                GetProtein = Lattice2(:,:,protein_num);
                GetProtein(:,1) = GetProtein(:,1)+p;
                for m = 1:4
                    GetProtein(m,:) = StayInGrid(GetProtein(m,:),SqDim);                  
                    if Filled(GetProtein(m,:),Lattice2) == 1
                        checkOccupancy = 1;
                    end    
                end    
                if checkOccupancy == 0%CheckOccupancy = 1 means selected location is occupied and not availaible to translate.
                    Lattice2(:,:,protein_num) = GetProtein(:,:);
                end    
            elseif directionXorY == 2%Selected protein will translate vertically
                    GetProtein = Lattice2(:,:,protein_num);
                    GetProtein(:,2) = GetProtein(:,2)+p;
                    for m = 1:4
                        GetProtein(m,:) = StayInGrid(GetProtein(m,:),SqDim);                  
                        if Filled(GetProtein(m,:),Lattice2) == 1
                            checkOccupancy = 1;
                        end    
                    end    
                    if checkOccupancy == 0%CheckOccupancy = 1 means selected location is occupied and not availaible to translate. 
                        Lattice2(:,:,protein_num) = GetProtein(:,:);
                    end           
            end    
        end
        if TypeOfMove == 1 || TypeOfMove == 2
            Index = randi(4,1,1);%randomly selecting from 4 residues  
            Pn = Lattice2(Index,:,protein_num);
            PsMv = [];%Possible Movement array
            if TypeOfMove == 1  %This condition means selected residue will do a corner move.
                if (Index == 1 || Index == 4)
                    continue;   %The end residues(1st and 4th of each protein) can't do corner move, so we skip to next iteration of protein.
                end     
                if Lattice2(Index-1,1,protein_num)==Lattice2(Index+1,1,protein_num) || Lattice2(Index-1,2,protein_num)==Lattice2(Index+1,2,protein_num) 
                    continue;   %If we select 2nd/3rdresidue and the x or y coordinate of both 1st and 3rd/2nd and 4th are same respectively, then we skip to next protein as here,it shows that the protein is linear around selected residue.
                end
                %If distance between selected and its preceding and successing residue are both same, then residue is making a corner move within the Lattice grid system
                c = Lattice2(Index-1,:,protein_num);  d = Lattice2(Index+1,:,protein_num);
                m = Lattice2(Index,:,protein_num);
                a = m(1)+(d(1)+c(1)-2*m(1));  b = m(2)+(d(2)+c(2)-2*m(2));

                if Filled([a b],Lattice2) == 0
                    PsMv(end+1,:) = [a b];
                end                            
            end       
            if TypeOfMove == 2   %This condition means selected residue will do a end move.
                if (Index == 2 || Index == 3)
                     continue;
                end    
                if (Index == 1)
                    p= 1;
                else 
                    p = -1;%means index = 4 
                end    
                if Pn(1) == Lattice2(Index + p,1,protein_num)%Here, both end residue and adjacent residue has same x coordinate in lattice.
                    a = Pn(1)-1; b = Pn(1)+1; c = Lattice2(Index + p,2,protein_num);
                    point = StayInGrid([a c],SqDim);point2 = StayInGrid([b c],SqDim); 
                    if Filled(point,Lattice2)==0
                        PsMv(end+1,:) = [a c];
                    end
                    if Filled(point2,Lattice2)==0
                        PsMv(end+1,:) = [b c];
                    end
                elseif Pn(2) == Lattice2(Index + p,2,protein_num)%Here, both end residue and adjacent residue has same y coordinate in lattice.
                    a = Pn(2)-1; b = Pn(2)+1;c = Lattice2(Index +p,1,protein_num);
                    point = StayInGrid([c a],SqDim);point2 = StayInGrid([c b],SqDim); 
                    if Filled(point,Lattice2)==0
                        PsMv(end+1,:) = [c a]; %#ok<*SAGROW>
                    end
                    if Filled(point2,Lattice2)==0
                        PsMv(end+1,:) = [c b];
                    end                   
                end
            end
         end
        if height(PsMv)~=0 %We make a random movement for the selected residue from all set of moves in PsMv if possible.
           ran = randi(height(PsMv),1,1);
           choosenMovement = PsMv(ran,:);
           choosenMovement = StayInGrid(choosenMovement,SqDim);
           Lattice2(Index,:,protein_num) = choosenMovement; %Assigning original to modified lattice if it is stable.
        end                  
    end             
    Ecalc = Energy_A(Lattice2,EInitial,SqDim,InteractionEnergy); %Calculates Lattice Energy 
    letsCompareEnergy = rand(1,1,1);
    BoltzmannWeight = exp(InitialSystemEnergy(end) - Ecalc);    
    if Ecalc < InitialSystemEnergy(end)%Applying Metropolis Monte Carlo Criteria for movement.
        Lattice(:,:,:) = Lattice2;
        InitialSystemEnergy(end+1) = Ecalc;
        StableEnergyIndex(end+1)=Iteration;
    elseif Ecalc >= InitialSystemEnergy(end)                
        if BoltzmannWeight > letsCompareEnergy
           Lattice(:,:,:) = Lattice2;
           InitialSystemEnergy(end+1) = Ecalc;
           StableEnergyIndex(end+1)= Iteration;
        end
    end
    if rem(Iteration,2000) == 0 %Plotting residues in lattice for each 2000 Iterations
        clf;hold on;        
        [X,Y] = meshgrid(0:1:SqDim,0:1:SqDim);
        plot(X,Y,'b.'); title("Number of Moves :"+Iteration);
        xlabel("Interaction Energy : "+InteractionEnergy+" Unit"+"         Lattice System Energy : "+InitialSystemEnergy(end)+" Units",'fontweight','bold','fontsize',12);
        grid on;grid minor;  ax = gca; ax.GridAlpha = 0.3; ax.MinorGridAlpha = 0.3;
        for i = 1:6
            if i == 1%We also check distance between residues between plotting, as during grid crossing , it should not show a covalent bond(by a connected line) in grid
                if distance(Lattice(1,:,i),Lattice(2,:,i)) == 1
                    hold on;plot(Lattice(1:2,1,i),Lattice(1:2,2,i),'b');                   
                end
                if distance(Lattice(2,:,i),Lattice(3,:,i)) == 1
                    hold on;plot(Lattice(2:3,1,i),Lattice(2:3,2,i),'b');                    
                end
                if distance(Lattice(3,:,i),Lattice(4,:,i)) == 1
                    hold on;plot(Lattice(3:4,1,i),Lattice(3:4,2,i),'b');                   
                end
                hold on;%Here, we plot each residue with specific color as per given in assignment. 
                h1 = plot(Lattice(1,1,i),Lattice(1,2,i),'-o',"MarkerFaceColor",'b','MarkerEdgeColor','b','DisplayName','A');hold on;          
                h2 = plot(Lattice(2,1,i),Lattice(2,2,i),'-o',"MarkerFaceColor",'g','MarkerEdgeColor','b','DisplayName','B');hold on;
                h3 = plot(Lattice(3,1,i),Lattice(3,2,i),'-o',"MarkerFaceColor",'y','MarkerEdgeColor','b','DisplayName','C');hold on;
                h4 = plot(Lattice(4,1,i),Lattice(4,2,i),'-o',"MarkerFaceColor",'r','MarkerEdgeColor','b','DisplayName','D');hold on;
                
            elseif i>1
                if distance(Lattice(1,:,i),Lattice(2,:,i)) == 1
                    hold on; plot(Lattice(1:2,1,i),Lattice(1:2,2,i),'b');                   
                end
                if distance(Lattice(2,:,i),Lattice(3,:,i)) == 1
                    hold on; plot(Lattice(2:3,1,i),Lattice(2:3,2,i),'b');                    
                end
                if distance(Lattice(3,:,i),Lattice(4,:,i)) == 1
                    hold on;  plot(Lattice(3:4,1,i),Lattice(3:4,2,i),'b');                   
                end
                hold on;
                plot(Lattice(1,1,i),Lattice(1,2,i),'-o',"MarkerFaceColor",'b','MarkerEdgeColor','b'); hold on;     
                plot(Lattice(2,1,i),Lattice(2,2,i),'-o',"MarkerFaceColor",'g','MarkerEdgeColor','b'); hold on;
                plot(Lattice(3,1,i),Lattice(3,2,i),'-o',"MarkerFaceColor",'y','MarkerEdgeColor','b'); hold on;
                plot(Lattice(4,1,i),Lattice(4,2,i),'-o',"MarkerFaceColor",'r','MarkerEdgeColor','b'); hold on;
            end
        end 
        drawnow;grid on;
    end
end
hl = legend( [ h1, h2, h3, h4 ], 'A', 'B', 'C', 'D','orientation', 'horizontal', 'location', 'Northoutside' );%It shows legends in figure 1 for residue/beads.
hold off;

figure(2);%This figure window shows variation of Lattice Energy with Number of Monte Carlo Iterations.
plot(StableEnergyIndex(:),InitialSystemEnergy(:),'-p')
xlabel('Monte Carlo Iterations','fontweight','bold','fontsize',12);
ylabel('Lattice Energy','fontweight','bold','fontsize',12);
a = width(InitialSystemEnergy)-1;
title("Number of Times Lattice Configuration Changed : " + a);
set(gcf, 'Position',  [860, 100, 600, 550]);
drawnow;

function isOccupied = Filled(E,LocalLattice)%This function check whether given point location is occupied.
    isOccupied = false;
    for j = 1:6
       for i = 1:4
          if E == LocalLattice(i,:,j)
              isOccupied = true;
          end
       end
    end        
end

function Energy1 = Energy_A(GivenLattice,EInitial,SqDim,InteractionEnergy)%GivenLattice refers to all the moves of present iteration applied previous Lattice
 Energy1 = EInitial;
 for i = 1:4
     for j = 1:6
         for k = j+1:6
             if distance(GivenLattice(i,:,j),GivenLattice(i,:,k)) == 1 || distance(GivenLattice(i,:,j),GivenLattice(i,:,k)) == SqDim + 1
                 if GivenLattice(i,1,j) - GivenLattice(i,1,k) == 0 || GivenLattice(i,2,j) - GivenLattice(i,2,k) == 0
                    Energy1 = Energy1 + InteractionEnergy;
                 end    
             end
         end     
     end    
 end   
end
%Energy2 calculates Energy of lattice system considering possible interactions between each pair of residues of different protein.
function Energy2 = Energy_B(GivenLattice,EInitial,SqDim,InteractionEnergy)
  Energy2 = EInitial;
  for j = 1:5
      A = [];
      for k = j+1:6  %error here
          for q=1:4
            A(end+1,:)= GivenLattice(q,:,k); %#ok<AGROW>
          end
      end
      for i = 1:4
         for m = 1:height(A)
             if distance(GivenLattice(i,:,j),A(m,:)) == 1 || distance(GivenLattice(i,:,j),A(m,:)) == 1 + SqDim
                 if GivenLattice(i,1,j) - A(m,1) == 0 || GivenLattice(i,2,j) - A(m,2) == 0
                    Energy2 = Energy2 + InteractionEnergy;
                 end
             end
         end    
      end    
  end    
end

%Energy3 calculates Energy of lattice system considering possible interaction between A - A And B - B only.
function Energy3 = Energy_C(GivenLattice,EInitial,SqDim,InteractionEnergy)
  Energy3 = EInitial;
 for i = 1:2
     for j = 1:6
         for k = j+1:6
             if distance(GivenLattice(i,:,j),GivenLattice(i,:,k)) == 1 || distance(GivenLattice(i,:,j),GivenLattice(i,:,k)) == SqDim + 1
                 if GivenLattice(i,1,j) - GivenLattice(i,1,k) == 0 || GivenLattice(i,2,j) - GivenLattice(i,2,k) == 0
                    Energy3 = Energy3 + InteractionEnergy;
                 end    
             end
         end     
     end    
 end   
end
%distance function calculates distance between any two point location on grid
function distance1 = distance(arr1,arr2)
    distance1 = ((arr1(1) - arr2(1)).^(2))+ ((arr1(2) - arr2(2)).^(2));
end    

function choosen = StayInGrid(choosenMovement,SqDim) 
    choosen = choosenMovement;
    if choosen(1)<0       
        choosen(1) = choosen(1) + 1 + SqDim;      
    elseif choosen(1)>SqDim
        choosen(1) = choosen(1) - 1 - SqDim;   
    elseif choosen(2)<0 
        choosen(2) = choosen(2) + 1 + SqDim;
    elseif choosen(2)>SqDim
        choosen(2) = choosen(2) - 1 - SqDim;
    end  
end   