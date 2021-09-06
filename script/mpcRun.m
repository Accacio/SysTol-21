clear
clc
simK=10;
Te=0.5;
n=5;

% options = optimset('Display', 'off')
opts = optimoptions('quadprog', 'MaxIter', 200, 'Display','off');

X0(:,1) = [18.3 3.0]';

umin = [0 ];
umax = [2 ];

Rf=[5];
Ri=[2.5];
Ro=[0.5];
Cres=[5];
Cs=[8];

Wt(:,:,1) = [20*ones(1,floor(simK/2)) 21*ones(1,simK-floor(simK/2))];

A=@(Cres,Cs,Rf,Ri,Ro) [-1/(Cres*Rf)-1/(Cres*Ri) 1/(Cres*Ri); 1/(Cs*Ri) -1/(Cs*Ro)-1/(Cs*Ri)];
B=@(Cres,Cs,Rf,Ro) [10/Cres 1/(Cres*Rf); 0 1/(Cs*Ro)];
Bp=@(Cres,Cs,Rf,Ro) [10/Cres; 0];
C=[1 0]; D=[0 0];

psys = c2d(ss(A(Cres,Cs,Rf,Ri,Ro),Bp(Cres,Cs,Rf,Ro),C,D(:,1)),Te);
sys.A = psys.A; sys.B = psys.B; sys.C = psys.C; sys.D = psys.D;
K=1;
xt(:,K) = X0;
Wt = Wt(:,:);
sys.umin = umin(:);  sys.umax = umax(:);

s= size(sys.A,1);   % number of states  of each subsystem
o= size(sys.C,1); % number of outputs of each subsystem
c= size(sys.B,2); % number of inputs  of each subsystem
Q(:,:)=eye(n*size(sys.C,1)); % no x no
R(:,:)=eye(n*size(sys.B,2)); % nc x nc

paren = @(x, varargin) x(varargin{:}); %
                                       % apply index in created elements
curly = @(x, varargin) x{varargin{:}}; %

getMmat=@(A,C,n) ...
 cell2mat( ...
    (repmat(mat2cell(C,1),1,n) .* ...
     (repmat(mat2cell(A,size(A,2)),1,n).^num2cell(1:n)) ...
    )'  ...
  );

getCmat=@(A, B, C, n) ...
    cell2mat( ...
      paren( ...
       (repmat(mat2cell(C, 1), 1, n+1) .* ...
         (horzcat(zeros(size(A)), repmat(mat2cell(A, size(A,2)),1,n).^num2cell(1:n))) .* ...
           repmat(mat2cell(B, size(B,1), size(B,2)), 1, n+1)) ...
           , tril(toeplitz(1:n))+1));


getH=@(Cmat,Q,R,n) round(Cmat'*Q*Cmat+R*eye(n),10);

getF=@(Cmat,Mmat,Q,xt,Wt) Cmat'*Q*(Mmat*xt-Wt);
getc=@(Mmat,Q,Xt,Wt) Xt'*Mmat'*Q*Mmat*Xt-2*Wt'*Q*Mmat*Xt+Wt'*Q*Wt;

Cmat(:,:) = getCmat(sys.A,sys.B,sys.C, n);
Mmat(:,:) = getMmat(sys.A,sys.C,n);

H(:,:) = getH(Cmat(:,:), Q, R, n);


for K=1:simK
    newWt(:,:)=kron(Wt(:,K),ones(1,n));
    F(:,:)=getF(Cmat(:,:), Mmat(:,:), Q(:,:), xt(:,K), newWt(:,:)');
    const(:,:)= getc(Mmat(:,:), Q,xt(:,K), newWt(:,:)');

    [u(:),J(:),~,~,l] = quadprog(H(:,:), F(:,:), [], [], [], [], umin(:)*ones(c*n,1), umax(:)*ones(c*n,1), [], opts);


    uHist(:,K)=u;
    xtPred(:,K)=Mmat*xt(:,K)+Cmat*u';
    xt(:,K+1) = sys.A(:,:)*xt(:,K)+sys.B(:,:)*u(1);
end

%%
for K=1:simK
    stairs(1:simK,Wt,'r')
    hold on
    plot(1:K,xt(1,1:K),'b*')
    stairs(K:K+n,[xt(1,K) xtPred(:,K)'],'k--')
    hold off
    xlim([1 simK])
    ylim([18 22])
    set(gcf, 'PaperPosition', [0 0 8 7])
    set(gca,'DefaultTextFontSize',18)
    set(gca,'FontSize',18)
    set(gca,'FontName','cmr18')
    set(gcf, 'PaperSize', [8 7])
    hx=get(gca,'xlabel');
    hy=get(gca,'ylabel');
    hz=get(gca,'zlabel');
    ht=get(gca,'title');
    hl=legend();
    lgd=legend({'Setpoint' 'State' 'Prediction'},'Location','southeast');
    lgd.FontSize = 16;
    lgd.FontWeight = 'bold';
    set(hx,'FontSize',18)
    set(hy,'FontSize',18)
    set(hz,'FontSize',18)
    set(hl,'FontSize',18)

    saveas(gcf,['../img/mpc' num2str(K) '.pdf'])
end
