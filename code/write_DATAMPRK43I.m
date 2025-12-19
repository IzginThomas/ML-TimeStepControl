function write_DATAMPRK43I()
clc, clear, close all
warning('off','all')
tic

solverstr = 'MPRK43Iadap';

%%%  Parameter of the method
solvpara1int = linspace(0.5,1,2); % 0.5 <= alpha <=1
solvpara2int = linspace(0.5,0.75,2);  % 0.5 <= beta <= 0.75
solvpara1int = [1];
solvpara2int = [0.5];
alpha0 = 1/6*(3 + (3-2*sqrt(2))^(1/3) + (3+2*sqrt(2))^(1/3)); % restricts alpha and beta, see below

%%% tolerances
tolint = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];
%tolint = [1e-4,1e-3,1e-2,1e-1];
%tolint = [1e-8];

%%% parameter of PID controller
dint = 0.01;
%dint = 0.2;
Beta1int = 0.1:dint:1;
Beta2int = -0.4:dint:0.2;
Beta3int = 0:dint:0.1;
%Beta3int = [0];
alpha2int = 1./(2:6);
kappa2 = 1;

%%% parameter of test problems PR3 and PR4
paraint = 0:0.1:0.5;
%paraint = [0.2];

%%% string of test problems
%tests = {{'PR3',paraint},{'PR4',paraint},'robertson','strat'}; 
tests ={'robertson','strat',{'PR4',paraint}};
%tests ={{'PR4',paraint}};

flagpool = true;
for testpair = tests
    testpair = testpair{:};
    if length(testpair) == 2
        flag = true;
        Test = testpair{1};
        testpara =testpair{2};
        lengthpara = length(testpara);
    else
        flag = false;
        Test = testpair;
        lengthpara = 1;
    end
    A = zeros(length(solvpara2int)*length(solvpara1int)*lengthpara*length(tolint)*length(Beta1int)*length(Beta2int)*length(Beta3int)* length(alpha2int),12);
    k = 1;
    for solvpara1 = solvpara1int
        for solvpara2 = solvpara2int
            if (1/3<=solvpara1 && solvpara1<2/3) && (2/3<= solvpara2 && solvpara2 <= 3*solvpara1*(1-solvpara1))
                paraflag = true;
            elseif (2/3<solvpara1 && solvpara1<alpha0) && (2/3>= solvpara2 && solvpara2 >= 3*solvpara1*(1-solvpara1))
                paraflag = true;
            elseif solvpara1 > alpha0 && (solvpara2 <= 2/3 && (3*solvpara1-2)/(6*solvpara1-3)<=solvpara2)
                paraflag = true;
            else
                paraflag = false;
            end
            if paraflag
                for beta1 = Beta1int
                    for beta2 = Beta2int
                        %if beta1 + beta2 > 0
                        for tol = tolint
                            for beta3 = Beta3int
                                for alpha2 = alpha2int
                                    if flag
                                        for para=testpara
                                            A(k,:) = [beta1, beta2, beta3, alpha2, tol, para,solvpara1,solvpara2, 0, 0, 0, 0];
                                            k = k+1;
                                        end
                                    else
                                        A(k,:) = [beta1, beta2, beta3, alpha2, tol, NaN ,solvpara1,solvpara2, 0, 0, 0, 0];
                                        k = k +1;
                                    end
                                end
                            end
                        end
                        %end
                    end
                end
            end
        end
    end


    A( all(A == 0,2), : ) = []; % removes zero rows from A
    lengthA = size(A,1);

    if flagpool
        %parpool(4);
        pool = parpool(60);
        parfor_opts = parforOptions(pool);
        flagpool = false;
    end

    q = parallel.pool.DataQueue;
    n_completed = 0;
    afterEach(q,@update_progress);
     parfor (i = 1:lengthA, parfor_opts)
    %parfor i = 1:lengthA
        %tic
        [Beta1, Beta2, Beta3, Alpha2, tol, para,solvpara1,solvpara2] = matsplit(A(i,:));
        [succstep, cnt_rej, err] = ODESolverTestSuite(Beta1, Beta2, Beta3, Alpha2, tol,Test,para,kappa2,solverstr,solvpara1,solvpara2);
        A(i,:) = [Beta1, Beta2, Beta3, Alpha2, tol, para,solvpara1,solvpara2, succstep, cnt_rej, err(1), err(2)];
        send(q,i);
        %toc
    end


    strfile = sprintf('DATA_%s_%s.txt',solverstr,Test);
    fileID = fopen(strfile,'w');

    %fprintf(fileID,'%12s %12s %12s %12s %12s %12s %12s\n','Beta1','Beta2', 'Beta3', 'tol','para','solvpara1','solvpara2' , 'accepted', 'rejected', 'error');
    fprintf(fileID,'%1.2f %1.2f %1.2f %1.2f %.e %1.1f %1.2f %1.2f %d %d %20.16f %20.16f\n',A');

    fclose(fileID);
end

%%% update function for displaying the progress in percent
function update_progress(~)
    n_completed = n_completed + 1;
    clc
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%s_%s:%4.2f\n','Progress MPRK43I',Test,n_completed*100/lengthA);
end
toc
end
