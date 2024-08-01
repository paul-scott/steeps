within SolarTherm.Validation.Models.PCM;

model Stefan_VariableMesh
    import Modelica.Math.Special.erf;
    import Modelica.Math.Special.erfc;
//Start of the HTF-PCM Model
    //Material Parameters
	parameter Real T_wall = 1118;//1078.0 "K";
	parameter Real T_init = 1028;//1068.0 "k";
    parameter Real L = 0.1 "m"; //Length of a plane wall with length L
    parameter Real A = 1.0 "m^2";
    parameter Integer n3 = 20; //20 //Total Number of nodes used in L for f3
    parameter Integer n2 = 30; //30 //Total Number of nodes used in L for f2
    parameter Integer n1 = 45; //45 //Total Number of nodes used in L for f1

    //Exact Solution Parameters;
    parameter Real pi = 2.0*Modelica.Math.asin(1.0);
    parameter Real St_s = 856.0*(1073.0-T_init)/4.312e5;
    parameter Real St_l = 1198.0*(T_wall-1073.0)/4.312e5;
    parameter Real v = sqrt((1.25/(2160.0*1198.0))/(4.0/(2160.0*856.0))); // defined at the root of ratio of thermal diffusivity of liquid/solid
    parameter Real alpha_l = 1.25/(2160.0*1198.0);
    parameter Real alpha_s = 4.0/(2160.0*856.0);
    parameter Real lamb = Lambda(St_s,St_l,v); //Solved via Newton-Raphson Method
    parameter Real r = 1.1 + 1.0e-6; //Mesh growth ratio, if you want 1.0, use 1.0 + 1.0e-6 to prevent div0 error
    Real X_exact (start = 0.0); //Location of melting front
    Real Q_exact (start = 0.0); //Total energy absorbed
    Real Slope_exact;
    

    //Discretization Variables f1
    parameter Real x1[n1] (start = Location(Thickness(r,L,n1),L,n1)) "m";//Distance along the x-axis, starts with incorrect value
    Real T1[n1] (start = fill(T_init,n1)) "K";//Temperature array
    Real h1[n1] (start = fill((T_init-300.0)*856.0,n1)); //Enthalpy
    Real k1[n1] (start = fill(4.0,n1)); //"W/mK"
    Real dx1[n1] (start=Thickness(r,L,n1)); //thickness of node
    Real Slope1; //Slope at the heating boundary
    parameter Real m1[n1] = Mass((Rho(0.0)*A)*A*L,Thickness(r,L,n1),n1); //mass of a single node
    Real f1[n1] (start = fill(0.0,n1));//liquid mass fraction
    Real X1 (start = 0.0); //Approx location of the melting front
    
    //Analytics

    Real T_exact1[n1]; //Exact temperature array    
    Real Q_tot1 (start = 0.0);

    Real Q_res1;
    Real T_res1[n1]; //Residual Temperature = Numerical - Exact
    Real Q_ErrorPct1; //Percentage error in heat absorbed
    Real X_res1;
    
    //Discretization Variables f2
    parameter Real x2[n2] (start = Location(Thickness(r,L,n2),L,n2)) "m";//Distance along the x-axis, starts with incorrect value
    Real T2[n2] (start = fill(T_init,n2)) "K";//Temperature array
    Real h2[n2] (start = fill((T_init-300.0)*856.0,n2)); //Enthalpy
    Real k2[n2] (start = fill(4.0,n2)); //"W/mK"
    Real dx2[n2] (start=Thickness(r,L,n2)); //thickness of node
    Real Slope2; //Slope at the heating boundary
    parameter Real m2[n2] = Mass(Rho(0.0)*A*L,Thickness(r,L,n2),n2); //mass of a single node
    Real f2[n2] (start = fill(0.0,n2));//liquid mass fraction
    Real X2 (start = 0.0); //Approx location of the melting front
    
    //Analytics

    Real T_exact2[n2]; //Exact temperature array    
    Real Q_tot2 (start = 0.0);

    Real Q_res2;
    Real T_res2[n2]; //Residual Temperature = Numerical - Exact
    Real Q_ErrorPct2; //Percentage error in heat absorbed
    Real X_res2;
    
    //Discretization Variables f3
    parameter Real x3[n3] (start = Location(Thickness(r,L,n3),L,n3)) "m";//Distance along the x-axis, starts with incorrect value
    Real T3[n3] (start = fill(T_init,n3)) "K";//Temperature array
    Real h3[n3] (start = fill((T_init-300.0)*856.0,n3)); //Enthalpy
    Real k3[n3] (start = fill(4.0,n3)); //"W/mK"
    Real dx3[n3] (start=Thickness(r,L,n3)); //thickness of node
    Real Slope3; //Slope at the heating boundary
    parameter Real m3[n3] = Mass(Rho(0.0)*A*L,Thickness(r,L,n3),n3); //mass of a single node
    Real f3[n3] (start = fill(0.0,n3));//liquid mass fraction
    Real X3 (start = 0.0); //Approx location of the melting front
    
    //Analytics

    Real T_exact3[n3]; //Exact temperature array    
    Real Q_tot3 (start = 0.0);

    Real Q_res3;
    Real T_res3[n3]; //Residual Temperature = Numerical - Exact
    Real Q_ErrorPct3; //Percentage error in heat absorbed
    Real X_res3;
    
    //Rate Of Convergence Calculations
    //Heat Absorbed
    Real ROC12_Q;
    Real ROC23_Q;
    Real ORC_Q;
    
    //Temperature Arrays
    Real ROC12_T;
    Real ROC23_T;
    Real ORC_T;
    
    //L2-norm of Temperature arrays
    Real Norm_T1;
    Real Norm_T2;
    Real Norm_T3;
    
    //Richardson extrapolation based on order of convergence
    //Temperature arrays are hard to process since they do not line up
    Real Q_Richardson;
    Real q_Richardson;
    Real Q_Richardson_ErrorPct; 
    
    //Slope Error
    Real Slope_Error1;
    Real Slope_Error2;
    Real Slope_Error3;
    Real Slope_ErrorPct1;
    Real Slope_ErrorPct2;
    Real Slope_ErrorPct3;   
    
    Real q_exact;
    Real q1;
    Real q2;
    Real q3;
    Real q_ErrorPct1;
    Real q_ErrorPct2;
    Real q_ErrorPct3;
    Real q_res1;
    Real q_res2;
    Real q_res3;
    Real ROC12_q;
    Real ROC23_q;
    Real ORC_q;
    
    
    Real Asymptote_T;
    Real Asymptote_Q;
    
    Real ORC_X;
    Real dx_min;
    Real ROC_X12;
    Real ROC_X23;
    
    Real M1;
    Real M2;
    Real M3;
    Real D1;
    Real D2;
    Real D3;
    
    Real X_pct1;
    Real X_pct2;
    Real X_pct3;
algorithm
    when max(max(T1[n1],T2[n2]),T3[n3]) > T_init + 0.01*(T_wall-T_init) then
      terminate("No longer semi-infinite");
    end when;
equation
//Independent of n
    dx_min = min(dx1);
    X_exact = 2.0*lamb*sqrt(alpha_l*time);
    Q_exact = (2160.0*4.312E5*St_l*X_exact)/(sqrt(pi)*lamb*erf(lamb)); //Multiplied by 2 due to the two HT sides
    q_exact = -1.0*Slope_exact*1.25*A;//Exact instantaneous heat transfer rate
    
    //Now for the exact temperature arrays
    //Set f1
    if time > 1.0 then
      Slope_exact = (1073.0-T_wall)/(sqrt(pi*alpha_l*time)*erf(lamb));
      for i in 1:n1 loop
        if x1[i] < X_exact then //we know it is liquid
          T_exact1[i] = T_wall - (T_wall-1073.0)*(erf(x1[i]/(2.0*sqrt(alpha_l*time)))/erf(lamb));
        elseif x1[i] > X_exact then //we know it is a solid
          T_exact1[i] = T_init + (1073.0-T_init)*(erfc(x1[i]/(2.0*sqrt(alpha_s*time)))/erfc(lamb*v));
        else //it is melting
          T_exact1[i] = 1073.0;
        end if;
      end for;
    
    //Set f2
      for i in 1:n2 loop
        if x2[i] < X_exact then //we know it is liquid
          T_exact2[i] = T_wall - (T_wall-1073.0)*(erf(x2[i]/(2.0*sqrt(alpha_l*time)))/erf(lamb));
        elseif x2[i] > X_exact then //we know it is a solid
          T_exact2[i] = T_init + (1073.0-T_init)*(erfc(x2[i]/(2.0*sqrt(alpha_s*time)))/erfc(lamb*v));
        else //it is melting
          T_exact2[i] = 1073.0;
        end if;
      end for;
    
    //Set f3
      for i in 1:n3 loop
        if x3[i] < X_exact then //we know it is liquid
          T_exact3[i] = T_wall - (T_wall-1073.0)*(erf(x3[i]/(2.0*sqrt(alpha_l*time)))/erf(lamb));
        elseif x3[i] > X_exact then //we know it is a solid
          T_exact3[i] = T_init + (1073.0-T_init)*(erfc(x3[i]/(2.0*sqrt(alpha_s*time)))/erfc(lamb*v));
        else //it is melting
          T_exact3[i] = 1073.0;
        end if;
      end for;
    else
      T_exact1 = fill(T_init,n1);
      T_exact2 = fill(T_init,n2);
      T_exact3 = fill(T_init,n3);
      Slope_exact = 0.0;
    end if;
    
    
    //Discretization f1
    for i in 1:n1 loop
      (T1[i], f1[i]) = Tf(h1[i]);
      k1[i] = K(f1[i]);
      dx1[i] = m1[i]/(Rho(f1[i])*A);
    end for;
    
    //Heat transfer in the middle section
	for i in 2:n1-1 loop
      m1[i]*der(h1[i])=A*((2.0*((k1[i+1]*k1[i])/(dx1[i]*k1[i+1]+dx1[i+1]*k1[i]))*(T1[i+1]-T1[i])) - (2.0*((k1[i]*k1[i-1])/(dx1[i]*k1[i-1]+dx1[i-1]*k1[i]))*(T1[i]-T1[i-1])));
    end for;
    
    //Left boundary i = 1
	m1[1]*der(h1[1])= -1.0*k1[1]*A*Slope1 - ((A*2.0*k1[1]*k1[2])/((k1[1]*dx1[2]+k1[2]*dx1[1])))*(T1[1]-T1[2]);
	
	//Right boundary i = n
	m1[n1]*der(h1[n1])= (2.0*((k1[n1]*k1[n1-1])/(dx1[n1]*k1[n1-1]+dx1[n1-1]*k1[n1]))*(T1[n1-1]-T1[n1]));
	
    Slope1 = (T1[1]-T_wall)/(0.5*dx1[1]);
	X1 = sum(f1.*dx1);
	der(Q_tot1) = -1.0*k1[1]*A*Slope1;
	
	//Exact Solution
    T_res1 = T1 - T_exact1;
    Q_res1 = Q_tot1 - Q_exact;
    X_res1 = X1 - X_exact;
    Q_ErrorPct1 = 100.0*(Q_tot1-Q_exact)/Q_exact;
    q1 = -1.0*k1[1]*A*Slope1;

	
    //Discretization f2
    for i in 1:n2 loop
      (T2[i], f2[i]) = Tf(h2[i]);
      k2[i] = K(f2[i]);
      dx2[i] = m2[i]/(Rho(f2[i])*A);
    end for;
    
    //Heat transfer in the middle section
	for i in 2:n2-1 loop
      m2[i]*der(h2[i])=A*((2.0*((k2[i+1]*k2[i])/(dx2[i]*k2[i+1]+dx2[i+1]*k2[i]))*(T2[i+1]-T2[i])) - (2.0*((k2[i]*k2[i-1])/(dx2[i]*k2[i-1]+dx2[i-1]*k2[i]))*(T2[i]-T2[i-1])));
    end for;
    
    //Left boundary i = 1
	m2[1]*der(h2[1])= -1.0*k2[1]*A*Slope2 - ((A*2.0*k2[1]*k2[2])/((k2[1]*dx2[2]+k2[2]*dx2[1])))*(T2[1]-T2[2]);
	
	//Right boundary i = n
	m2[n2]*der(h2[n2])= (2.0*((k2[n2]*k2[n2-1])/(dx2[n2]*k2[n2-1]+dx2[n2-1]*k2[n2]))*(T2[n2-1]-T2[n2]));
	
	Slope2 = (T2[1]-T_wall)/(0.5*dx2[1]);
	X2 = sum(f2.*dx2);
	der(Q_tot2) = -1.0*k2[1]*A*Slope2;
	
	//Exact Solution
    T_res2 = T2 - T_exact2;
    Q_res2 = Q_tot2 - Q_exact;
    X_res2 = X2 - X_exact;
    Q_ErrorPct2 = 100.0*(Q_tot2-Q_exact)/Q_exact;
    q2 = -1.0*k2[1]*A*Slope2;

	
	//Discretization f2
    for i in 1:n3 loop
      (T3[i], f3[i]) = Tf(h3[i]);
      k3[i] = K(f3[i]);
      dx3[i] = m3[i]/(Rho(f3[i])*A);
    end for;
    
    //Heat transfer in the middle section
	for i in 2:n3-1 loop
      m3[i]*der(h3[i])=A*((2.0*((k3[i+1]*k3[i])/(dx3[i]*k3[i+1]+dx3[i+1]*k3[i]))*(T3[i+1]-T3[i])) - (2.0*((k3[i]*k3[i-1])/(dx3[i]*k3[i-1]+dx3[i-1]*k3[i]))*(T3[i]-T3[i-1])));
    end for;
    
    //Left boundary i = 1
	m3[1]*der(h3[1])= -1.0*k3[1]*A*Slope3 - ((A*2.0*k3[1]*k3[2])/((k3[1]*dx3[2]+k3[2]*dx3[1])))*(T3[1]-T3[2]);
	
	//Right boundary i = n
	m3[n3]*der(h3[n3])= (2.0*((k3[n3]*k3[n3-1])/(dx3[n3]*k3[n3-1]+dx3[n3-1]*k3[n3]))*(T3[n3-1]-T3[n3]));
	
	Slope3 = (T3[1]-T_wall)/(0.5*dx3[1]);
	X3 = sum(f3.*dx3);
	der(Q_tot3) = -1.0*k3[1]*A*Slope3;
	
	//Exact Solution
    T_res3 = T3 - T_exact3;
    Q_res3 = Q_tot3 - Q_exact;
    X_res3 = X3 - X_exact;
    Q_ErrorPct3 = 100.0*(Q_tot3-Q_exact)/Q_exact;
    q3 = -1.0*k3[1]*A*Slope3;
    q_res1 = q1 - q_exact;
    q_res2 = q2 - q_exact;
    q_res3 = q3 - q_exact;

    
    //Rate of Convergence
	if time > 1.0 then
      ROC12_Q = log(abs(Q_res1/Q_res2))/log(n2/n1); //because step size is inversely related to n
      ROC23_Q = log(abs(Q_res2/Q_res3))/log(n3/n2);
      ROC12_q = log(abs(q_res1/q_res2))/log(n2/n1);
      ROC23_q = log(abs(q_res2/q_res3))/log(n3/n2);
      ORC_q = log(abs((q3-q2)/(q2-q1)))/log(n1/n2);
      if abs(Q_tot2-Q_tot1) > 100.0 then
        ORC_Q = log(abs((Q_tot3-Q_tot2)/(Q_tot2-Q_tot1)))/log(n1/n2); //Assuming n1/n2=n2/n3=r
        Q_Richardson = Q_tot1 + (Q_tot1 - Q_tot2)/((n1/n2)^ORC_Q-1.0);
        Q_Richardson_ErrorPct = 100.0*abs(Q_Richardson - Q_exact)/Q_exact;
      else
        ORC_Q = 0.0;
        Q_Richardson = 0.0;
        Q_Richardson_ErrorPct = 0.0;
      end if;
      q_Richardson = q1 + (q1 - q2)/((n1/n2)^ORC_q-1.0);
      if abs((X2-X1)) == 0.0 then
        ORC_X = 0.0;
      else
        if abs(X3-X2) == 0.0 then
          ORC_X = 0.0;
        else
          ORC_X = log(abs((X3-X2)/(X2-X1)))/log(n1/n2);
        end if;
      end if;
      //Need the L2-norm formula is sum of percentage errors divided by number of points
      Norm_T1 = sqrt(sum((T_res1.*T_res1)./(T_exact1.*T_exact1)))/n1;
      Norm_T2 = sqrt(sum((T_res2.*T_res2)./(T_exact2.*T_exact2)))/n2;
      Norm_T3 = sqrt(sum((T_res3.*T_res3)./(T_exact3.*T_exact3)))/n3;
      
      ROC12_T = log(abs(Norm_T1/Norm_T2))/log(n2/n1);
      ROC23_T = log(abs(Norm_T2/Norm_T3))/log(n3/n2);
      if (Norm_T2-Norm_T1) > 1E-12 then
        ORC_T = log(abs((Norm_T3-Norm_T2)/(Norm_T2-Norm_T1)))/log(n1/n2);
      else
        ORC_T = 0.0;
      end if;
      
      q_ErrorPct1 = 100.0*(q1-q_exact)/q_exact;
      q_ErrorPct2 = 100.0*(q2-q_exact)/q_exact;
      q_ErrorPct3 = 100.0*(q3-q_exact)/q_exact;
      //Slopes
      //Slope Error
      Slope_Error1 = Slope1 - Slope_exact;
      Slope_Error2 = Slope2 - Slope_exact;
      Slope_Error3 = Slope3 - Slope_exact;
      Slope_ErrorPct1 = 100.0*((Slope1 - Slope_exact)/Slope_exact);
      Slope_ErrorPct2 = 100.0*((Slope2 - Slope_exact)/Slope_exact);
      Slope_ErrorPct3 = 100.0*((Slope3 - Slope_exact)/Slope_exact);
      
      Asymptote_T= ((Norm_T2-Norm_T3)*Norm_T1)/(Norm_T2*(Norm_T1-Norm_T2));
      Asymptote_Q= ((Q_tot2-Q_tot3)*Q_tot1)/(Q_tot2*(Q_tot1-Q_tot2));
      ROC_X23 = log(abs(X_res2/X_res3))/log(n3/n2);
      ROC_X12 = log(abs(X_res1/X_res2))/log(n2/n1);
      
    else
      ROC12_Q = 0;
      ROC23_Q = 0;
      ORC_Q = 0;
      
      ROC12_q = 0;
      ROC23_q = 0;
      ORC_q = 0;
      
      Norm_T1 = sqrt(sum((T_res1.*T_res1)./(T_exact1.*T_exact1)))/n1;
      Norm_T2 = sqrt(sum((T_res2.*T_res2)./(T_exact2.*T_exact2)))/n2;
      Norm_T3 = sqrt(sum((T_res3.*T_res3)./(T_exact3.*T_exact3)))/n3;
      
      ROC12_T = 0;
      ROC23_T = 0;
      ORC_T = 0;
      
      Q_Richardson = 0;
      Q_Richardson_ErrorPct = 0;
      q_Richardson = 0;
      
      Slope_Error1 = 0;
      Slope_Error2 = 0;
      Slope_Error3 = 0;
      Slope_ErrorPct1 = 0;
      Slope_ErrorPct2 = 0;
      Slope_ErrorPct3 = 0;
      
      Asymptote_T= 0.0;
      Asymptote_Q= 0.0;
      
      q_ErrorPct1 = 0.0;
      q_ErrorPct2 = 0.0;
      q_ErrorPct3 = 0.0;
      ORC_X = 0;
      ROC_X12 = 0;
      ROC_X23 = 0;
    end if;
    
    M1 = sum(m1);
    M2 = sum(m2);
    M3 = sum(m3);
    D1 = sum(dx1);
    D2 = sum(dx2);
    D3 = sum(dx3);
    
    X_pct1 = 100.0*X_res1/X_exact;
    X_pct2 = 100.0*X_res2/X_exact;
    X_pct3 = 100.0*X_res3/X_exact;
    
    annotation(experiment(StopTime = 330, StartTime = 0, Tolerance = 1e-6, Interval = 0.1));
end Stefan_VariableMesh;