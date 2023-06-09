%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <cgmeehan@alumni.caltech.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

        function [mlY, mlScrPos]=javaPointToMatLab(javaX, javaY, screens)
            if nargin<3
                screens=javaScreens;
            end
            N=length(screens);            
            for i=1:N
                pe=screens{i};
                if pe.y==0 && pe.x==0
                    defaultIdx=i;
                    defaultPe=pe;
                end
            end
            idx=0;
            for i=1:N
                pe=screens{i};
                if javaX>=pe.x &&  javaX<pe.x+pe.width
                    if javaY>=pe.y && javaY<pe.y+pe.height
                        idx=i;
                        break;
                    end
                end
            end
            if idx==0
                mlScrPos=[];
                 mlY=javaY;
                 return;
            end
            pe=screens{idx};
            mlTopY=defaultPe.height-pe.y;
            tmp=(defaultPe.height-javaY)-mlTopY;
            mlY=mlTopY+tmp;
            idx=0;
            if idx==defaultIdx
                mlScrPos=[1, 1, pe.width, pe.height];
            else
                if pe.y<0
                    y=defaultPe.height-(pe.height+pe.y);
                else
                    y=(defaultPe.height-pe.y)-pe.height;
                end
                mlScrPos=[pe.x+1, y+1, pe.width, pe.height];
            end
            
        end
