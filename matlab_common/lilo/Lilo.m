function [A,Xmozaic]=stitching(mozaic, files)
%% Stitching of N tiles
% Written by Maxime Deforet, April 5. 2012.
% INPUT :
% - mozaic : array of "overlapping graph"
% exemple : [[1,2];[2,3];[3,4];[4,5];[3,6];[6,7];[7,8];[8,9];[9,10]]; 
% for 
%      2  1
%    3  4  5
%    6  7  8
%     10  9
% - files : array of N cells, each of these cells is the full path of the tile i.
%     example : 
%     folder='/home/your_images';
%     files=dir(folder)
%     files=sort_nat({files.name});
%     for i=1:10
%         Files{i}=fullfile(folder,files{i});
%     end
% 
% OUTPUT : 
% - A : stitched image from tiles
% - Xmozaic : xy positions of tiles

a=double(imread(files{1}));

%% définition des Dx et Dy. 
Ntiles=length(files);
% "path of stitching" :
% mozaic=[[1,2];[2,3];[3,4];[4,5];[3,6];[6,7];[7,8];[8,9];[9,10]];


Xmozaic=zeros(Ntiles,2);
DXmozaic=zeros(Ntiles-1,2);
BiCor=xcorrf2(ones(size(a)),ones(size(a)))/(size(a,1)*size(a,2)); % normalisation des effets de taille finie
BW=zeros(size(BiCor));
BW(:,1:15)=1;BW(:,end-15:end)=1;BW(1:15,:)=1;BW(end-15:end,:)=1; %supprime violemment les pics aux bords.
for i=1:size(mozaic,1)
    a=double(imread(files{mozaic(i,1)}));
    b=double(imread(files{mozaic(i,2)}));
    mag=imtophat(xcorrf2(a,b)./BiCor,ones(5)); %imtophat pour enlever les correlations a grande portée, et laisser le tout petit pics
    mag(BW==1)=0;
    [fX fY] = find((mag == (max(max(mag)))));
    DXmozaic(i,:)=[fX fY];
end
DXmozaic(:,1)=DXmozaic(:,1)-size(a,1);
DXmozaic(:,2)=DXmozaic(:,2)-size(a,2);
for i=2:Ntiles
    Xmozaic(i,:)=Xmozaic(mozaic(find(mozaic(:,2)==i,1,'first'),1),:)+DXmozaic(find(mozaic(:,2)==i,1,'first'),:);
end


%% stitching : reconstruction de l'image totale avec les coordonnées du stitching obtenues ci-dessus.    
[sx,sy]=size(a);
[ax,ay]=meshgrid(1:sy,1:sx);

xx=[];
yy=[];
aa=[];
for i=1:Ntiles
    a=double(imread(files{i}));
    xx=[xx;ax(:)+Xmozaic(i,2)];
    yy=[yy;ay(:)+Xmozaic(i,1)];
    aa=[aa;a(:)]; %liste des valeurs des pixels
end
xx=xx-min(xx)+1;
yy=yy-min(yy)+1;
ii=sub2ind([max(yy),max(xx)],yy,xx); %liste des pixels des images stitchées

[iii,idx]=sort(ii);
aa_sorted=aa(idx);
% iii et aa_sorted sont les nouveaux indices et pixels, ok.

fin=circshift(logical((circshift(iii,-1)-iii)),1)-logical((circshift(iii,-1)-iii)); %fin des zones de recouvrement
deb=circshift(logical((circshift(iii,1)-iii)),-1)-logical((circshift(iii,1)-iii)); %debut des zones de recouvrement dans la liste des indices
[occ,occ2]=rude(iii); %occ est le nombre d'occurence de chaque indice / occ2 est la valeur de chacun de ces indices
occi=rude(occ,occ); % occi, de la meme taille que iii, donne le nombre d'occurence en face des indices, memes doublés.
cum=cumsum(aa_sorted./rude(occ,occ)'); % somme cumulative, coeeficientée par le nombre d'occurence
%cum(fin==-1)-cum(find(deb==-1)-1)
[B,M,N]=unique(iii);
idx_unique=(M(occi(M)==1)); % : indice des uniques, dans l'espace de la liste des uniques (numéro de positions)
idx_superpos=iii(M(occi(M)~=1)); % indice des superpositions dans l'espace reel, lineaire.
iiii=[iii(idx_unique);idx_superpos];
aaaa=[aa_sorted(idx_unique);cum(fin==-1)-cum(find(deb==-1)-1)];

% On regroupe tout dans une sparse :
[xxx,yyy]=ind2sub([max(yy),max(xx)],iiii);
s=sparse(xxx(:),yyy(:),aaaa);

% Matrice reelle :
A=full(s);
Xmozaic=fliplr(Xmozaic);


function c = xcorrf2(a,b,pad)
%  c = xcorrf2(a,b)
%   Two-dimensional cross-correlation using Fourier transforms.
%       XCORRF2(A,B) computes the crosscorrelation of matrices A and B.
%       XCORRF2(A) is the autocorrelation function.
%       This routine is functionally equivalent to xcorr2 but usually faster.
%       See also XCORR2.

%       Author(s): R. Johnson
%       $Revision: 1.0 $  $Date: 1995/11/27 $

  if nargin<3
    pad='yes';
  end
  
  
  [ma,na] = size(a);
  if nargin == 1
    %       for autocorrelation
    b = a;
  end
  [mb,nb] = size(b);
  %       make reverse conjugate of one array
  b = conj(b(mb:-1:1,nb:-1:1));
  
  if strcmp(pad,'yes');
    %       use power of 2 transform lengths
    mf = 2^nextpow2(ma+mb);
    nf = 2^nextpow2(na+nb);
    at = fft2(b,mf,nf);
    bt = fft2(a,mf,nf);
  elseif strcmp(pad,'no');
    at = fft2(b);
    bt = fft2(a);
  else
    disp('Wrong input to XCORRF2'); return
  end
  
  %       multiply transforms then inverse transform
  c = ifft2(at.*bt);
  %       make real output for real input
  if ~any(any(imag(a))) && ~any(any(imag(b)))
    c = real(c);
  end
  %  trim to standard size
  
  if strcmp(pad,'yes');
    c(ma+mb:mf,:) = [];
    c(:,na+nb:nf) = [];
  elseif strcmp(pad,'no');
    c=fftshift(c(1:end-1,1:end-1));
    
%    c(ma+mb:mf,:) = [];
%    c(:,na+nb:nf) = [];
  end
  
  
  
  
  % VEC       = rude(LEN,VAL)
%		run-length DEcoding
%
% [LEN,VAL] = rude(VEC)
%		run-length ENcoding
%
% P         = rude;
%		retrieve subroutine handles in structure P
%   P.d		for run-length DEcoding
%			VEC = P.d(LEN,VAL)
%   P.e		for run-length ENcoding
%			[LEN,VAL] = P.e(VEC)
%
% LEN	: repeat each VAL corresponding LEN times to create VEC
% VAL	: 1xN template array of VALs to be repeated LEN times
%	  - numericals
%	  - strings
%	  - cells (any contents)
% VEC	: DEcode = reconstruced output vector from LEN/VAL
%	  ENcode = input vector to be encoded into LEN/VAL
%
% NOTE
% 1:	LEN <= 0 will remove corresponding VALs
% 2:	<NaN>s and <+Inf>s are treated as not-equal (expected behavior!)
%
% USAGE EXAMPLE
%	vec = rude([1 2 3 0 4],[10 inf nan pi 20])
% %	vec = 10 Inf Inf NaN NaN NaN 20 20 20 20
%	[len,val] = rude(vec)
% %	len = 1 1 1 1 1 1 4 % note nan~=nan / inf~=inf!
% %	val = 10 Inf Inf NaN NaN NaN 20
%
%	s = rude;
%	v.x = pi;
%	w.x = pi; % note <v> and <w> are equal!
%	vec = s.d([1 0 3 2 2 2 3],{'a' 'b' 'cd' 1:3 v w magic(3)})
% %	vec = 'a' 'cd' 'cd' 'cd' 1x3D 1x3D 1x1S 1x1S 1x1S 1x1S 3x3D
%	[len,val] = s.e(vec)
% %	len = 1 3 2 4 3
% %	val = 'a' 'cd' 1x3D 1x1S 3x3D

% created:
%	us	18-Nov-2004
% modified:
%	us	30-Nov-2004 20:42:03	/ TMW FEX

function	[p1,p2]=rude(varargin)

		p2=[];
	if	~nargin & nargout
		p1.d=@rl_decode;
		p1.e=@rl_encode;
	elseif	~nargin
		help(mfilename);
		return;
	else
	if	nargin == 1
		[p1,p2]=rl_encode(varargin{1});
	elseif	nargin >= 2
		p1=rl_decode(varargin{1:2});
	end
	end
		return;
%--------------------------------------------------------------------------------
% run-length decoder
function	vec=rl_decode(len,val)

		lx=len>0 & ~(len==inf);
	if	~any(lx)
		vec=[];
		return;
	end
	if	numel(len) ~= numel(val)
		error(...
		sprintf(['rl-decoder: length mismatch\n',...
			 'len = %-1d\n',...
			 'val = %-1d'],...
			  numel(len),numel(val)));
	end
		len=len(lx);
		val=val(lx);
		val=val(:).';
		len=len(:);
		lc=cumsum(len);
		lx=zeros(1,lc(end));
		lx([1;lc(1:end-1)+1])=1;
		lc=cumsum(lx);
		vec=val(lc);
		return;
%--------------------------------------------------------------------------------
% run-length encoder
function	[len,val]=rl_encode(vec)

	switch	class(vec)
	case	'cell'
		[len,val]=rl_encodec(vec);
	case	'char'
		[len,val]=rl_encoden(double(vec));
		val=char(val);
	otherwise
		[len,val]=rl_encoden(vec);
	end
		return;
%--------------------------------------------------------------------------------
% run-length encode doubles
function	[len,val]=rl_encoden(vec)

	if	isempty(vec)
		len=0;
		val=[];
		return;
	end
		vec=vec(:).';
		vx=[1 diff(double(vec))];
		vx=vx~=0;
		val=vec(vx);
		vc=cumsum(vx).';
		len=accumarray(vc,ones(size(vc))).';
		return;
%--------------------------------------------------------------------------------
% run-length encode cells
function	[len,val]=rl_encodec(vec)

		len=0;
		val={};
		vl=length(vec)+1;
		cl=cellfun('length',vec);
		tix=0;
		tlen=1;
	for	i=1:length(cl)
	if	cl(i)
		tmpl=vec{i};
	for	j=i+1:vl
	if	j>length(cl) || ~isequalwithequalnans(tmpl,vec{j})
		tix=tix+1;
		val{tix}=tmpl;
		len(tix)=tlen;
		cl(i:j-1)=0;
		tlen=1;
		break;
	else
		tlen=tlen+1;
	end
	end
	end
	end
		return;
%--------------------------------------------------------------------------
%------
