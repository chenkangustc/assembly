	  �-  �   k820309    ?          14.0        iZ                                                                                                           
       E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\debugassembly\driving_output.f90 DRIVING_OUTPUT          RUN_OUTPUT                                                 
                     @               �                '           #FRIC    #GEOM    #MESH    #PROPERTY $   #BOUNDARY 2   #INITDATA @   #CONFACTOR_ Q   #POW Z   #THERMAL ^   #ALLOC i   #CLEAN l   #SET p   #INIT u            �                                           	            �                                             #ASSMGEOM                  @                               '            #RFUEL    #GASGAP    #SHELLTHICK    #ASSMSHELLTHICK 	   #ACROSSFLAT 
   #HEIGHT    #N_PIN    #SET             �                                           	            �                                          	            �                                          	            �                               	           	            �                               
           	            �                                          	            �                                             1     �   � $                      �                    #SET_ASSMGEOM    #     @     @                                                #THIS    #RFUEL    #GASGAP    #SHELLTHICK    #ASSMSHELLTHICK    #ACROSSFLAT    #HEIGHT    #N_PIN          
                                            #ASSMGEOM          
                                       	        
                                       	        
                                       	        
                                       	        
                                       	        
                                       	        
                                                  �                                              #ASSMMESH                  @                               '            #NF    #NG    #NS    #NY    #SET             �                                                       �                                                      �                                                      �                                             1     �   � $                      �                    #SET_ASSMMESH    #     @     @                                                #THIS    #NF     #NG !   #NS "   #NY #         
                                            #ASSMMESH          
                                                
                                  !             
                                  "             
                                  #                �                               $     �   0      #MATERIAL %                 @              @           %     '�            #RHO &   #SHC '   #CTC (   #DVS )   #HTC *   #INIT +          �                               &                  	        &           &                              �                               '        0         	        &           &                              �                               (        `         	        &           &                              �                               )        �         	        &           &                              �                               *        �         	        &                       1     �   � $                      �      +              #INIT_MATERIAL ,   #     @     @                            ,                    #THIS -   #NF .   #NG /   #NS 0   #NY 1         
                                -     �       #MATERIAL %         
                                  .             
                                  /             
                                  0             
                                  1                �                               2             #BOUNDARY 3                 @                          3     '            #TIN 4   #TOUT 5   #UIN 6   #UOUT 7   #PIN 8   #POUT 9   #INIT :            �                               4            	            �                               5           	            �                               6           	            �                               7           	            �                               8           	            �                               9           	   1     �   � $                      �      :              #INIT_BOUNDARY ;   #     @     @                            ;                    #THIS <   #TIN =   #UIN >   #PIN ?         
                                <            #BOUNDARY 3         
                                  =     	        
                                  >     	        
                                  ?     	           �                               @        ,     #ASSMINIT A                 @                          A     '            #TI B   #PI C   #UI D   #TIN E   #PIN F   #UIN G   #SET H            �                               B            	            �                               C           	            �                               D           	            �                               E           	            �                               F           	            �                               G           	   1     �   � $                      �      H              #SET_ASSMINIT I   #     @     @                            I                    #THIS J   #TI K   #PI L   #UI M   #TIN N   #PIN O   #UIN P         
                                J            #ASSMINIT A         
                                  K     	        
                                  L     	        
                                  M     	        
                                  N     	        
                                  O     	        
                                  P     	           �                               Q        D     #CONFACTOR R                 @                          R     '            #ALPHA S   #SIGMA T   #SET U            �                               S            	            �                               T           	   1     �   � $                      �      U              #SET_CONFACTOR V   #     @     @                            V                    #THIS W   #ALPHA X   #SIGMA Y         
                                W            #CONFACTOR R         
                                  X     	        
                                  Y     	           �                               Z     H   L     #ASSMPOW [                 @              @           [     'H            #POWER \   #FQ_CORE ]          �                               \                  	        &                              �                               ]        $         	        &                                �                               ^     x   �  	   #THERMAL _                 @              @           _     'x            #TEMPERATURE `   #PRESSURE a   #VELOCITY b   #INIT c          �                               `                  	        &           &                              �                               a        0         	        &                              �                               b        T         	        &                       1     �   � $                      �      c              #INIT_THERMAL d   #     @     @                            d                    #THIS e   #TEMPERATURE f   #PRESSURE g   #VELOCITY h         
                                e     x       #THERMAL _         
                                  f     	        
                                  g     	        
                                  h     	  1     �   � $                      �      i         
     #ALLOC_ASSEMBLY j   #     @     @                            j                    #THIS k         
                                k           #SYS_ASSEMBLY    1     �   � $                      �      l              #FREE_ASSEMBLY m   #     @     @                            m                   #FREE_ASSEMBLY%ALLOCATED n   #THIS o                                           n     ALLOCATED       
                                o           #SYS_ASSEMBLY    1     �   � $                      �      p              #SET_ASSEMBLY q   #     @     @                            q                    #THIS r   #REINPUTDATA s         
                                r           #SYS_ASSEMBLY          
                                  s     X      #SYS_RE_INPUT t   1     �   � $                      �      u              #INIT_ASSEMBLY v   #     @     @                            v                    #THIS w         
                                w           #SYS_ASSEMBLY    #     @                                   x                                   @                           t     'X            #NF y   #NG z   #NS {   #NY |   #NPIN }   #XF ~   #XG    #XS �   #XOS �   #ACF �   #HEIGHT �   #F �   #POUT �   #FLOWIN �   #TIN �   #UIN �   #PIN �   #TI �   #UI �   #PI �   #ALPHA �   #SIGMA �   #SET �            �                               y                        �                               z                       �                               {                       �                               |                       �                               }                       �                               ~           	            �                                          	            �                               �           	            �                               �         	   	            �                               �     $   
   	            �                               �     (      	            �                               �     ,      	            �                               �     0      	            �                               �     4      	            �                               �     8      	            �                               �     <      	            �                               �     @      	            �                               �     D      	            �                               �     H      	            �                               �     L      	            �                               �     P      	            �                               �     T      	   1     �   � $                      �      �              #SET_INPUTDATA �   #     @     @                            �                    #THIS �                                         �     X       #SYS_RE_INPUT t      �   �      fn#fn $   4     b   uapp(DRIVING_OUTPUT    K  <   J  ASSM_GLOBAL 1   �  �       SYS_ASSEMBLY+SYS_ASSEMBLY_HEADER 6   b  @   a   SYS_ASSEMBLY%FRIC+SYS_ASSEMBLY_HEADER 6   �  N   a   SYS_ASSEMBLY%GEOM+SYS_ASSEMBLY_HEADER )   �  �       ASSMGEOM+SYS_ASSM_HEADER /   �  @   a   ASSMGEOM%RFUEL+SYS_ASSM_HEADER 0   �  @   a   ASSMGEOM%GASGAP+SYS_ASSM_HEADER 4     @   a   ASSMGEOM%SHELLTHICK+SYS_ASSM_HEADER 8   _  @   a   ASSMGEOM%ASSMSHELLTHICK+SYS_ASSM_HEADER 4   �  @   a   ASSMGEOM%ACROSSFLAT+SYS_ASSM_HEADER 0   �  @   a   ASSMGEOM%HEIGHT+SYS_ASSM_HEADER /     @   a   ASSMGEOM%N_PIN+SYS_ASSM_HEADER -   _  R   a   ASSMGEOM%SET+SYS_ASSM_HEADER -   �  �      SET_ASSMGEOM+SYS_ASSM_HEADER 2   a  J   a   SET_ASSMGEOM%THIS+SYS_ASSM_HEADER 3   �  8   a   SET_ASSMGEOM%RFUEL+SYS_ASSM_HEADER 4   �  8   a   SET_ASSMGEOM%GASGAP+SYS_ASSM_HEADER 8     8   a   SET_ASSMGEOM%SHELLTHICK+SYS_ASSM_HEADER <   S  8   a   SET_ASSMGEOM%ASSMSHELLTHICK+SYS_ASSM_HEADER 8   �  8   a   SET_ASSMGEOM%ACROSSFLAT+SYS_ASSM_HEADER 4   �  8   a   SET_ASSMGEOM%HEIGHT+SYS_ASSM_HEADER 3   �  8   a   SET_ASSMGEOM%N_PIN+SYS_ASSM_HEADER 6   3  N   a   SYS_ASSEMBLY%MESH+SYS_ASSEMBLY_HEADER )   �  m       ASSMMESH+SYS_ASSM_HEADER ,   �  @   a   ASSMMESH%NF+SYS_ASSM_HEADER ,   .	  @   a   ASSMMESH%NG+SYS_ASSM_HEADER ,   n	  @   a   ASSMMESH%NS+SYS_ASSM_HEADER ,   �	  @   a   ASSMMESH%NY+SYS_ASSM_HEADER -   �	  R   a   ASSMMESH%SET+SYS_ASSM_HEADER -   @
  n      SET_ASSMMESH+SYS_ASSM_HEADER 2   �
  J   a   SET_ASSMMESH%THIS+SYS_ASSM_HEADER 0   �
  8   a   SET_ASSMMESH%NF+SYS_ASSM_HEADER 0   0  8   a   SET_ASSMMESH%NG+SYS_ASSM_HEADER 0   h  8   a   SET_ASSMMESH%NS+SYS_ASSM_HEADER 0   �  8   a   SET_ASSMMESH%NY+SYS_ASSM_HEADER :   �  N   a   SYS_ASSEMBLY%PROPERTY+SYS_ASSEMBLY_HEADER )   &  {       MATERIAL+SYS_ASSM_HEADER -   �  |   a   MATERIAL%RHO+SYS_ASSM_HEADER -     |   a   MATERIAL%SHC+SYS_ASSM_HEADER -   �  |   a   MATERIAL%CTC+SYS_ASSM_HEADER -     |   a   MATERIAL%DVS+SYS_ASSM_HEADER -   �  l   a   MATERIAL%HTC+SYS_ASSM_HEADER .   �  S   a   MATERIAL%INIT+SYS_ASSM_HEADER .   P  n      INIT_MATERIAL+SYS_ASSM_HEADER 3   �  J   a   INIT_MATERIAL%THIS+SYS_ASSM_HEADER 1     8   a   INIT_MATERIAL%NF+SYS_ASSM_HEADER 1   @  8   a   INIT_MATERIAL%NG+SYS_ASSM_HEADER 1   x  8   a   INIT_MATERIAL%NS+SYS_ASSM_HEADER 1   �  8   a   INIT_MATERIAL%NY+SYS_ASSM_HEADER :   �  N   a   SYS_ASSEMBLY%BOUNDARY+SYS_ASSEMBLY_HEADER )   6  �       BOUNDARY+SYS_ASSM_HEADER -   �  @   a   BOUNDARY%TIN+SYS_ASSM_HEADER .   �  @   a   BOUNDARY%TOUT+SYS_ASSM_HEADER -   =  @   a   BOUNDARY%UIN+SYS_ASSM_HEADER .   }  @   a   BOUNDARY%UOUT+SYS_ASSM_HEADER -   �  @   a   BOUNDARY%PIN+SYS_ASSM_HEADER .   �  @   a   BOUNDARY%POUT+SYS_ASSM_HEADER .   =  S   a   BOUNDARY%INIT+SYS_ASSM_HEADER .   �  i      INIT_BOUNDARY+SYS_ASSM_HEADER 3   �  J   a   INIT_BOUNDARY%THIS+SYS_ASSM_HEADER 2   C  8   a   INIT_BOUNDARY%TIN+SYS_ASSM_HEADER 2   {  8   a   INIT_BOUNDARY%UIN+SYS_ASSM_HEADER 2   �  8   a   INIT_BOUNDARY%PIN+SYS_ASSM_HEADER :   �  N   a   SYS_ASSEMBLY%INITDATA+SYS_ASSEMBLY_HEADER )   9  �       ASSMINIT+SYS_ASSM_HEADER ,   �  @   a   ASSMINIT%TI+SYS_ASSM_HEADER ,   �  @   a   ASSMINIT%PI+SYS_ASSM_HEADER ,   9  @   a   ASSMINIT%UI+SYS_ASSM_HEADER -   y  @   a   ASSMINIT%TIN+SYS_ASSM_HEADER -   �  @   a   ASSMINIT%PIN+SYS_ASSM_HEADER -   �  @   a   ASSMINIT%UIN+SYS_ASSM_HEADER -   9  R   a   ASSMINIT%SET+SYS_ASSM_HEADER -   �  �      SET_ASSMINIT+SYS_ASSM_HEADER 2     J   a   SET_ASSMINIT%THIS+SYS_ASSM_HEADER 0   V  8   a   SET_ASSMINIT%TI+SYS_ASSM_HEADER 0   �  8   a   SET_ASSMINIT%PI+SYS_ASSM_HEADER 0   �  8   a   SET_ASSMINIT%UI+SYS_ASSM_HEADER 1   �  8   a   SET_ASSMINIT%TIN+SYS_ASSM_HEADER 1   6  8   a   SET_ASSMINIT%PIN+SYS_ASSM_HEADER 1   n  8   a   SET_ASSMINIT%UIN+SYS_ASSM_HEADER <   �  O   a   SYS_ASSEMBLY%CONFACTOR_+SYS_ASSEMBLY_HEADER *   �  c       CONFACTOR+SYS_ASSM_HEADER 0   X  @   a   CONFACTOR%ALPHA+SYS_ASSM_HEADER 0   �  @   a   CONFACTOR%SIGMA+SYS_ASSM_HEADER .   �  S   a   CONFACTOR%SET+SYS_ASSM_HEADER .   +  d      SET_CONFACTOR+SYS_ASSM_HEADER 3   �  K   a   SET_CONFACTOR%THIS+SYS_ASSM_HEADER 4   �  8   a   SET_CONFACTOR%ALPHA+SYS_ASSM_HEADER 4     8   a   SET_CONFACTOR%SIGMA+SYS_ASSM_HEADER 5   J  M   a   SYS_ASSEMBLY%POW+SYS_ASSEMBLY_HEADER (   �  \       ASSMPOW+SYS_ASSM_HEADER .   �  l   a   ASSMPOW%POWER+SYS_ASSM_HEADER 0   _  l   a   ASSMPOW%FQ_CORE+SYS_ASSM_HEADER 9   �  M   a   SYS_ASSEMBLY%THERMAL+SYS_ASSEMBLY_HEADER (     {       THERMAL+SYS_ASSM_HEADER 4   �  |   a   THERMAL%TEMPERATURE+SYS_ASSM_HEADER 1     l   a   THERMAL%PRESSURE+SYS_ASSM_HEADER 1   {  l   a   THERMAL%VELOCITY+SYS_ASSM_HEADER -   �  R   a   THERMAL%INIT+SYS_ASSM_HEADER -   9   {      INIT_THERMAL+SYS_ASSM_HEADER 2   �   I   a   INIT_THERMAL%THIS+SYS_ASSM_HEADER 9   �   8   a   INIT_THERMAL%TEMPERATURE+SYS_ASSM_HEADER 6   5!  8   a   INIT_THERMAL%PRESSURE+SYS_ASSM_HEADER 6   m!  8   a   INIT_THERMAL%VELOCITY+SYS_ASSM_HEADER 7   �!  T   a   SYS_ASSEMBLY%ALLOC+SYS_ASSEMBLY_HEADER 3   �!  N      ALLOC_ASSEMBLY+SYS_ASSEMBLY_HEADER 8   G"  N   a   ALLOC_ASSEMBLY%THIS+SYS_ASSEMBLY_HEADER 7   �"  S   a   SYS_ASSEMBLY%CLEAN+SYS_ASSEMBLY_HEADER 2   �"  k      FREE_ASSEMBLY+SYS_ASSEMBLY_HEADER <   S#  >      FREE_ASSEMBLY%ALLOCATED+SYS_ASSEMBLY_HEADER 7   �#  N   a   FREE_ASSEMBLY%THIS+SYS_ASSEMBLY_HEADER 5   �#  R   a   SYS_ASSEMBLY%SET+SYS_ASSEMBLY_HEADER 1   1$  _      SET_ASSEMBLY+SYS_ASSEMBLY_HEADER 6   �$  N   a   SET_ASSEMBLY%THIS+SYS_ASSEMBLY_HEADER =   �$  N   a   SET_ASSEMBLY%REINPUTDATA+SYS_ASSEMBLY_HEADER 6   ,%  S   a   SYS_ASSEMBLY%INIT+SYS_ASSEMBLY_HEADER 2   %  N      INIT_ASSEMBLY+SYS_ASSEMBLY_HEADER 7   �%  N   a   INIT_ASSEMBLY%THIS+SYS_ASSEMBLY_HEADER    &  D       RUN_OUTPUT 1   _&        SYS_RE_INPUT+SYS_RE_INPUT_HEADER 4   r'  @   a   SYS_RE_INPUT%NF+SYS_RE_INPUT_HEADER 4   �'  @   a   SYS_RE_INPUT%NG+SYS_RE_INPUT_HEADER 4   �'  @   a   SYS_RE_INPUT%NS+SYS_RE_INPUT_HEADER 4   2(  @   a   SYS_RE_INPUT%NY+SYS_RE_INPUT_HEADER 6   r(  @   a   SYS_RE_INPUT%NPIN+SYS_RE_INPUT_HEADER 4   �(  @   a   SYS_RE_INPUT%XF+SYS_RE_INPUT_HEADER 4   �(  @   a   SYS_RE_INPUT%XG+SYS_RE_INPUT_HEADER 4   2)  @   a   SYS_RE_INPUT%XS+SYS_RE_INPUT_HEADER 5   r)  @   a   SYS_RE_INPUT%XOS+SYS_RE_INPUT_HEADER 5   �)  @   a   SYS_RE_INPUT%ACF+SYS_RE_INPUT_HEADER 8   �)  @   a   SYS_RE_INPUT%HEIGHT+SYS_RE_INPUT_HEADER 3   2*  @   a   SYS_RE_INPUT%F+SYS_RE_INPUT_HEADER 6   r*  @   a   SYS_RE_INPUT%POUT+SYS_RE_INPUT_HEADER 8   �*  @   a   SYS_RE_INPUT%FLOWIN+SYS_RE_INPUT_HEADER 5   �*  @   a   SYS_RE_INPUT%TIN+SYS_RE_INPUT_HEADER 5   2+  @   a   SYS_RE_INPUT%UIN+SYS_RE_INPUT_HEADER 5   r+  @   a   SYS_RE_INPUT%PIN+SYS_RE_INPUT_HEADER 4   �+  @   a   SYS_RE_INPUT%TI+SYS_RE_INPUT_HEADER 4   �+  @   a   SYS_RE_INPUT%UI+SYS_RE_INPUT_HEADER 4   2,  @   a   SYS_RE_INPUT%PI+SYS_RE_INPUT_HEADER 7   r,  @   a   SYS_RE_INPUT%ALPHA+SYS_RE_INPUT_HEADER 7   �,  @   a   SYS_RE_INPUT%SIGMA+SYS_RE_INPUT_HEADER 5   �,  S   a   SYS_RE_INPUT%SET+SYS_RE_INPUT_HEADER 2   E-  N      SET_INPUTDATA+SYS_RE_INPUT_HEADER 7   �-  N   a   SET_INPUTDATA%THIS+SYS_RE_INPUT_HEADER 