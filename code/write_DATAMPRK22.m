function write_DATAMPRK22()
clc, clear, close all
warning('off','all')
tic

solverstr = 'MPRK22adap';

%%%  Parameter of the method
%solvpara1int = [1/2,1]; %alpha
solvpara1int = [1];

%%% tolerances
tolint = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];
%tolint = [1e-4,1e-3,1e-2,1e-1];
%tolint = [1e-4];

%%% parameter of PID controller
dint = 0.01;
%dint = 0.1;
Beta1int = 0.1:dint:1;
Beta2int = -0.4:dint:0.2;
Beta3int = 0:dint:0.1;
%Beta3int = [0];

alpha2int = 1./(2:6);
%alpha2int = [1/6];
kappa2 = 1;

%%% parameter of test problems PR3 and PR4
paraint = 0:0.1:0.5;
%paraint = [0.2];

%%% string of test problems
%tests = {{'PR3',paraint},{'PR4',paraint},'robertson','strat'}; 
tests ={'robertson','strat',{'PR4',paraint}};
%tests ={'robertson',{'PR3',paraint}};
%tests = {{'PR3',paraint}};
    
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
    A = zeros(length(solvpara1int)*lengthpara*length(tolint)*length(Beta1int)*length(Beta2int)*length(Beta3int)* length(alpha2int),11);
    k = 1;
    for solvpara1 = solvpara1int
        for beta1 = Beta1int
            for beta2 = Beta2int
                % if or(beta1 + beta2 <= 0, beta2>-0.05)
                for tol = tolint
                    for beta3 = Beta3int
                        for alpha2 = alpha2int
                            if flag
                                for para=testpara
                                    A(k,:) = [beta1, beta2, beta3, alpha2, tol, para,solvpara1, 0, 0, 0, 0];
                                    k = k+1;
                                end
                            else
                                A(k,:) = [beta1, beta2, beta3, alpha2, tol, NaN ,solvpara1, 0, 0, 0, 0];
                                k = k +1;
                            end
                        end
                    end
                end
                %end
            end
        end
    end
    A( all(A == 0,2), : ) = []; % removes zero rows from A

    lengthA = size(A,1);
    %toc
    
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
  % parfor i = 1:lengthA
       %disp(i)
        %tic
        [Beta1, Beta2, Beta3, Alpha2, tol, para,solvpara1] = matsplit(A(i,:));
        [succstep, cnt_rej, err] = ODESolverTestSuite(Beta1, Beta2, Beta3, Alpha2, tol,Test,para,kappa2,solverstr,solvpara1); %err=[L2_rel_err_sol, L2_rel_err_secondinvariant]
        A(i,:) = [Beta1, Beta2, Beta3, Alpha2, tol, para, solvpara1, succstep, cnt_rej, err(1), err(2)];
        send(q,i);
        %toc
    end



    strfile = sprintf('DATA_%s_%s.txt',solverstr,Test);
    fileID = fopen(strfile,'w');

    %fprintf(fileID,'%12s %12s %12s %12s %12s %12s %12s\n','Beta1','Beta2', 'Beta3', 'tol','para','solvpara1' , 'accepted', 'rejected', 'error_sol', 'error_inv');
    fprintf(fileID,'%1.2f %1.2f %1.2f %1.2f %.e %1.1f %1.3f %d %d %20.16f %20.16f\n',A');

    fclose(fileID);
end



%%% update function for displaying the progress in percent
function update_progress(~)
    n_completed = n_completed + 1;
    clc
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%s_%s:%4.2f\n','Progress MPRK22',Test,n_completed*100/lengthA);
end
toc
end