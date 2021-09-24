% function P=drawZ2D(xc,R,C, t, z_flag)
%
% Affiche un polygone convexe repr�sentant la projection du zonotope Z(R) dans le plan 
% form� par les vecteurs [1 0 0...0]' et [0 1 0...0]' et translat�e de xc
%
% Le zonotope Z(R) repr�sente l'image du domaine [-1;1]^p par l'application lin�aire
% correspondant � la matrice R (n lignes et p colonnes):
%    Soit R=[R1 ... Rp] o� les Rj (j=1...p) sont des vecteurs colonne de d �l�ments
%    Z(R)={s1*R1+...+sp*Rp, |sj|<1, j=1...p}
% Remarque: Z(R) correspond � un domaine centr� sur z�ro
%
% xc (matrice d*1) : vecteur colonne d�signant le centre du zonotope � afficher
%                    Si xc=0 (scalaire), le zonotope affich� est centr� sur l'origine
% C                : couleur (exemple: 'b' pour bleu)
%                    si C='line', seul le contour est trac�
%                    si C='line,r', seul le contour est trac� en rouge
% t (scalaire)     : Si non vide, trac� de la projection 2D (x1,x2) dans un rep�re 3D (x,y,z)
%                    avec x=t, y=x1, z=x2. Utiliser rotate3D ou view pour modifier l'orientation. 
% z_flag           : Si t est non vide et z_flag=1, trac� de la projection 2D (x1,x2) dans un
%                    rep�re 3D (x,y,z) avec x=x1, y=x2, z=t. 
% P (matrice m*2)  : Liste des points d�finissant l'enveloppe (convexe) de la projection
%                    du zonotope Z(R) dans le plan P. Les points, tri�s par angle croissant,
%                    sont rep�r�s par leur abscisse (colonne 1) et leur ordonn�e (colonne 2).
%                    La liste de points correspond � un domaine centr� sur l'origine (z�ro).

% drawZ2D : from plotZ2D_V06_, C. Combastel, 2017-02

function P=drawZ2D(xc,R,C, t, z_flag)

% Initialisations et gestion des param�tres
% -----------------------------------------
if nargin<=2, C='b';     end   % Couleur      par d�faut
if nargin<=3, t=[];      end   % Instant/cote par d�faut
if nargin<=4, z_flag=0;  end   % Rotation 3D  par d�faut

n=size(R,1);    % Nombre de lignes de R
p=size(R,2);    % Nombre de colonnes de R (i.e. nombre de segments g�n�rateurs)
if xc==0
   xc=zeros(n,1);
end   
R  = R(1:2,:);  % Projection du zonotope dans le plan engendr� par [1 0 0...0]' et [0 1 0...0]'
xc = xc(1:2);   % Projection du centre

% Calcul 
% ------
R10=(R(1,:)==0);
R=R.*(ones(2,1)*(sign(R(1,:))+R10)); 
[AR,I]=sort(R(2,:)./(R(1,:)+eps*R10)); R=R(:,I);
p=size(R,2);
s=ones(p,1); P=R*s;
for i=1:(p-1),
   s(i)=-1;
   P=[P R*s];
end
P=[P -P]';

% Affichage du polygone correspondant � la projection 2D du zonotope
% ------------------------------------------------------------------
X=xc(1)+P(:,1);
Y=xc(2)+P(:,2);
if isempty(t)
   % Trac� de la projection 2D dans un rep�re 2D
   switch length(C),
      case 1, fill(X,Y,C);
      case 4, line([X;X(1)],[Y;Y(1)]);
      case 6, set(line([X;X(1)],[Y;Y(1)]),'color',C(6));
      otherwise, error('Format de C inconnu dans ''plotZ2D''');
   end
else
   % Trac� de la projection 2D dans un rep�re 3D
   Z=ones(size(P,1),1)*t;
   if (z_flag == 1)
      % Trac� de la projection 2D (x1,x2) dans un rep�re 3D (x,y,z) avec x=x1, y=x2, z=t
      xlabel('x_1');   ylabel('x_2');   zlabel('t');
      view([50 20]);   grid on;
   else   
      % Trac� de la projection 2D (x1,x2) dans un rep�re 3D (x,y,z) avec x=t,  y=x1, z=x2
      X_ = X;    Y_ = Y;    Z_ = Z;
      X  = Z_;   Y  = X_;   Z  = Y_;   
      xlabel('t');   ylabel('x_1');   zlabel('x_2');
      view([50 20]);   grid on;
   end
   switch length(C),
      case 1, fill3(X,Y,Z,C);
      case 4, line([X;X(1)],[Y;Y(1)],[Z;Z(1)]);
      case 6, set(line([X;X(1)],[Y;Y(1)],[Z;Z(1)]),'color',C(6));
      otherwise, error('Format de C inconnu dans ''plotZ2D''');
   end   
end  

return

