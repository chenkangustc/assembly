	  �-  �   k820309    ?          14.0        �Z                                                                                                           
       E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\debugassembly\assm_global.f90 ASSM_GLOBAL                  @                              
                     @               �                '           #FRIC    #GEOM    #MESH    #PROPERTY $   #TH_BOUNDARY 2   #INITDATA @   #CONFACTOR_ Q   #POW Z   #THERMAL ^   #ALLOC i   #CLEAN l   #SET p   #INIT u            �                                           	            �                                             #ASSMGEOM                  @                               '            #RFUEL    #GASGAP    #SHELLTHICK    #ASSMSHELLTHICK 	   #ACROSSFLAT 
   #HEIGHT    #N_PIN    #SET             �                                           	            �                                          	            �                                          	            �                               	           	            �                               
           	            �                                          	            �                                             1     �   � $                      �                    #SET_ASSMGEOM    #     @     @                                                #THIS    #RFUEL    #GASGAP    #SHELLTHICK    #ASSMSHELLTHICK    #ACROSSFLAT    #HEIGHT    #N_PIN          
                                            #ASSMGEOM          
                                       	        
                                       	        
                                       	        
                                       	        
                                       	        
                                       	        
                                                  �                                              #ASSMMESH                  @                               '            #NF    #NG    #NS    #NY    #SET             �                                                       �                                                      �                                                      �                                             1     �   � $                      �                    #SET_ASSMMESH    #     @     @                                                #THIS    #NF     #NG !   #NS "   #NY #         
                                            #ASSMMESH          
                                                
                                  !             
                                  "             
                                  #                �                               $     �   0      #MATERIAL %                 @              @           %     '�            #RHO &   #SHC '   #CTC (   #DVS )   #HTC *   #INIT +          �                               &                  	        &           &                              �                               '        0         	        &           &                              �                               (        `         	        &           &                              �                               )        �         	        &           &                              �                               *        �         	        &                       1     �   � $                      �      +              #INIT_MATERIAL ,   #     @     @                            ,                    #THIS -   #NF .   #NG /   #NS 0   #NY 1         
                                -     �       #MATERIAL %         
                                  .             
                                  /             
                                  0             
                                  1                �                               2             #TH_BOUNDARY 3                 @                          3     '            #P 4   #U 8   #T 9   #INIT :            �                               4               #BOUNDARY 5                 @                          5     '            #INLET 6   #OUTLET 7            �                               6            	            �                               7           	            �                               8              #BOUNDARY 5            �                               9              #BOUNDARY 5   1     �   � $                      �      :              #INIT_TH_BOUNDARY ;   #     @     @                            ;                    #THIS <   #TIN =   #UIN >   #PIN ?         
                                <            #TH_BOUNDARY 3         
                                  =     	        
                                  >     	        
                                  ?     	           �                               @        ,     #ASSMINIT A                 @                          A     '            #TI B   #PI C   #UI D   #TIN E   #PIN F   #UIN G   #SET H            �                               B            	            �                               C           	            �                               D           	            �                               E           	            �                               F           	            �                               G           	   1     �   � $                      �      H              #SET_ASSMINIT I   #     @     @                            I                    #THIS J   #TI K   #PI L   #UI M   #TIN N   #PIN O   #UIN P         
                                J            #ASSMINIT A         
                                  K     	        
                                  L     	        
                                  M     	        
                                  N     	        
                                  O     	        
                                  P     	           �                               Q        D     #CONFACTOR R                 @                          R     '            #ALPHA S   #SIGMA T   #SET U            �                               S            	            �                               T           	   1     �   � $                      �      U              #SET_CONFACTOR V   #     @     @                            V                    #THIS W   #ALPHA X   #SIGMA Y         
                                W            #CONFACTOR R         
                                  X     	        
                                  Y     	           �                               Z     H   L     #ASSMPOW [                 @              @           [     'H            #POWER \   #FQ_CORE ]          �                               \                  	        &                              �                               ]        $         	        &                                �                               ^     x   �  	   #THERMAL _                 @              @           _     'x            #TEMPERATURE `   #PRESSURE a   #VELOCITY b   #INIT c          �                               `                  	        &           &                              �                               a        0         	        &                              �                               b        T         	        &                       1     �   � $                      �      c              #INIT_THERMAL d   #     @     @                            d                    #THIS e   #TEMPERATURE f   #PRESSURE g   #VELOCITY h         
                                e     x       #THERMAL _         
                                  f     	        
                                  g     	        
                                  h     	  1     �   � $                      �      i         
     #ALLOC_ASSEMBLY j   #     @     @                            j                    #THIS k         
                                k           #SYS_ASSEMBLY    1     �   � $                      �      l              #FREE_ASSEMBLY m   #     @     @                            m                   #FREE_ASSEMBLY%ALLOCATED n   #THIS o                                           n     ALLOCATED       
                                o           #SYS_ASSEMBLY    1     �   � $                      �      p              #SET_ASSEMBLY q   #     @     @                            q                    #THIS r   #REINPUTDATA s         
                                r           #SYS_ASSEMBLY          
                                  s     X      #SYS_RE_INPUT t   1     �   � $                      �      u              #INIT_ASSEMBLY v   #     @     @                            v                    #THIS w         
                                w           #SYS_ASSEMBLY    #     @     @                            x                    #THIS y                                         y     X       #SYS_RE_INPUT t                                            z       #SYS_ASSEMBLY                  @                           t     'X            #NF {   #NG |   #NS }   #NY ~   #NPIN    #XF �   #XG �   #XS �   #XOS �   #ACF �   #HEIGHT �   #F �   #POUT �   #FLOWIN �   #TIN �   #UIN �   #PIN �   #TI �   #UI �   #PI �   #ALPHA �   #SIGMA �   #SET �            �                               {                        �                               |                       �                               }                       �                               ~                       �                                                      �                               �           	            �                               �           	            �                               �           	            �                               �         	   	            �                               �     $   
   	            �                               �     (      	            �                               �     ,      	            �                               �     0      	            �                               �     4      	            �                               �     8      	            �                               �     <      	            �                               �     @      	            �                               �     D      	            �                               �     H      	            �                               �     L      	            �                               �     P      	            �                               �     T      	   1     �   � $                      �      �              #SET_INPUTDATA x      �   �      fn#fn $   .  <   J   SYS_ASSEMBLY_HEADER 1   j  �       SYS_ASSEMBLY+SYS_ASSEMBLY_HEADER 6   H  @   a   SYS_ASSEMBLY%FRIC+SYS_ASSEMBLY_HEADER 6   �  N   a   SYS_ASSEMBLY%GEOM+SYS_ASSEMBLY_HEADER )   �  �       ASSMGEOM+SYS_ASSM_HEADER /   �  @   a   ASSMGEOM%RFUEL+SYS_ASSM_HEADER 0   �  @   a   ASSMGEOM%GASGAP+SYS_ASSM_HEADER 4     @   a   ASSMGEOM%SHELLTHICK+SYS_ASSM_HEADER 8   E  @   a   ASSMGEOM%ASSMSHELLTHICK+SYS_ASSM_HEADER 4   �  @   a   ASSMGEOM%ACROSSFLAT+SYS_ASSM_HEADER 0   �  @   a   ASSMGEOM%HEIGHT+SYS_ASSM_HEADER /     @   a   ASSMGEOM%N_PIN+SYS_ASSM_HEADER -   E  R   a   ASSMGEOM%SET+SYS_ASSM_HEADER -   �  �      SET_ASSMGEOM+SYS_ASSM_HEADER 2   G  J   a   SET_ASSMGEOM%THIS+SYS_ASSM_HEADER 3   �  8   a   SET_ASSMGEOM%RFUEL+SYS_ASSM_HEADER 4   �  8   a   SET_ASSMGEOM%GASGAP+SYS_ASSM_HEADER 8     8   a   SET_ASSMGEOM%SHELLTHICK+SYS_ASSM_HEADER <   9  8   a   SET_ASSMGEOM%ASSMSHELLTHICK+SYS_ASSM_HEADER 8   q  8   a   SET_ASSMGEOM%ACROSSFLAT+SYS_ASSM_HEADER 4   �  8   a   SET_ASSMGEOM%HEIGHT+SYS_ASSM_HEADER 3   �  8   a   SET_ASSMGEOM%N_PIN+SYS_ASSM_HEADER 6     N   a   SYS_ASSEMBLY%MESH+SYS_ASSEMBLY_HEADER )   g  m       ASSMMESH+SYS_ASSM_HEADER ,   �  @   a   ASSMMESH%NF+SYS_ASSM_HEADER ,   	  @   a   ASSMMESH%NG+SYS_ASSM_HEADER ,   T	  @   a   ASSMMESH%NS+SYS_ASSM_HEADER ,   �	  @   a   ASSMMESH%NY+SYS_ASSM_HEADER -   �	  R   a   ASSMMESH%SET+SYS_ASSM_HEADER -   &
  n      SET_ASSMMESH+SYS_ASSM_HEADER 2   �
  J   a   SET_ASSMMESH%THIS+SYS_ASSM_HEADER 0   �
  8   a   SET_ASSMMESH%NF+SYS_ASSM_HEADER 0     8   a   SET_ASSMMESH%NG+SYS_ASSM_HEADER 0   N  8   a   SET_ASSMMESH%NS+SYS_ASSM_HEADER 0   �  8   a   SET_ASSMMESH%NY+SYS_ASSM_HEADER :   �  N   a   SYS_ASSEMBLY%PROPERTY+SYS_ASSEMBLY_HEADER )     {       MATERIAL+SYS_ASSM_HEADER -   �  |   a   MATERIAL%RHO+SYS_ASSM_HEADER -     |   a   MATERIAL%SHC+SYS_ASSM_HEADER -     |   a   MATERIAL%CTC+SYS_ASSM_HEADER -   �  |   a   MATERIAL%DVS+SYS_ASSM_HEADER -   w  l   a   MATERIAL%HTC+SYS_ASSM_HEADER .   �  S   a   MATERIAL%INIT+SYS_ASSM_HEADER .   6  n      INIT_MATERIAL+SYS_ASSM_HEADER 3   �  J   a   INIT_MATERIAL%THIS+SYS_ASSM_HEADER 1   �  8   a   INIT_MATERIAL%NF+SYS_ASSM_HEADER 1   &  8   a   INIT_MATERIAL%NG+SYS_ASSM_HEADER 1   ^  8   a   INIT_MATERIAL%NS+SYS_ASSM_HEADER 1   �  8   a   INIT_MATERIAL%NY+SYS_ASSM_HEADER =   �  Q   a   SYS_ASSEMBLY%TH_BOUNDARY+SYS_ASSEMBLY_HEADER ,     c       TH_BOUNDARY+SYS_ASSM_HEADER .   �  N   a   TH_BOUNDARY%P+SYS_ASSM_HEADER )   �  [       BOUNDARY+SYS_ASSM_HEADER /   +  @   a   BOUNDARY%INLET+SYS_ASSM_HEADER 0   k  @   a   BOUNDARY%OUTLET+SYS_ASSM_HEADER .   �  N   a   TH_BOUNDARY%U+SYS_ASSM_HEADER .   �  N   a   TH_BOUNDARY%T+SYS_ASSM_HEADER 1   G  V   a   TH_BOUNDARY%INIT+SYS_ASSM_HEADER 1   �  i      INIT_TH_BOUNDARY+SYS_ASSM_HEADER 6     M   a   INIT_TH_BOUNDARY%THIS+SYS_ASSM_HEADER 5   S  8   a   INIT_TH_BOUNDARY%TIN+SYS_ASSM_HEADER 5   �  8   a   INIT_TH_BOUNDARY%UIN+SYS_ASSM_HEADER 5   �  8   a   INIT_TH_BOUNDARY%PIN+SYS_ASSM_HEADER :   �  N   a   SYS_ASSEMBLY%INITDATA+SYS_ASSEMBLY_HEADER )   I  �       ASSMINIT+SYS_ASSM_HEADER ,   �  @   a   ASSMINIT%TI+SYS_ASSM_HEADER ,   	  @   a   ASSMINIT%PI+SYS_ASSM_HEADER ,   I  @   a   ASSMINIT%UI+SYS_ASSM_HEADER -   �  @   a   ASSMINIT%TIN+SYS_ASSM_HEADER -   �  @   a   ASSMINIT%PIN+SYS_ASSM_HEADER -   	  @   a   ASSMINIT%UIN+SYS_ASSM_HEADER -   I  R   a   ASSMINIT%SET+SYS_ASSM_HEADER -   �  �      SET_ASSMINIT+SYS_ASSM_HEADER 2     J   a   SET_ASSMINIT%THIS+SYS_ASSM_HEADER 0   f  8   a   SET_ASSMINIT%TI+SYS_ASSM_HEADER 0   �  8   a   SET_ASSMINIT%PI+SYS_ASSM_HEADER 0   �  8   a   SET_ASSMINIT%UI+SYS_ASSM_HEADER 1     8   a   SET_ASSMINIT%TIN+SYS_ASSM_HEADER 1   F  8   a   SET_ASSMINIT%PIN+SYS_ASSM_HEADER 1   ~  8   a   SET_ASSMINIT%UIN+SYS_ASSM_HEADER <   �  O   a   SYS_ASSEMBLY%CONFACTOR_+SYS_ASSEMBLY_HEADER *     c       CONFACTOR+SYS_ASSM_HEADER 0   h  @   a   CONFACTOR%ALPHA+SYS_ASSM_HEADER 0   �  @   a   CONFACTOR%SIGMA+SYS_ASSM_HEADER .   �  S   a   CONFACTOR%SET+SYS_ASSM_HEADER .   ;  d      SET_CONFACTOR+SYS_ASSM_HEADER 3   �  K   a   SET_CONFACTOR%THIS+SYS_ASSM_HEADER 4   �  8   a   SET_CONFACTOR%ALPHA+SYS_ASSM_HEADER 4   "  8   a   SET_CONFACTOR%SIGMA+SYS_ASSM_HEADER 5   Z  M   a   SYS_ASSEMBLY%POW+SYS_ASSEMBLY_HEADER (   �  \       ASSMPOW+SYS_ASSM_HEADER .     l   a   ASSMPOW%POWER+SYS_ASSM_HEADER 0   o  l   a   ASSMPOW%FQ_CORE+SYS_ASSM_HEADER 9   �  M   a   SYS_ASSEMBLY%THERMAL+SYS_ASSEMBLY_HEADER (   (  {       THERMAL+SYS_ASSM_HEADER 4   �  |   a   THERMAL%TEMPERATURE+SYS_ASSM_HEADER 1     l   a   THERMAL%PRESSURE+SYS_ASSM_HEADER 1   �  l   a   THERMAL%VELOCITY+SYS_ASSM_HEADER -   �  R   a   THERMAL%INIT+SYS_ASSM_HEADER -   I   {      INIT_THERMAL+SYS_ASSM_HEADER 2   �   I   a   INIT_THERMAL%THIS+SYS_ASSM_HEADER 9   !  8   a   INIT_THERMAL%TEMPERATURE+SYS_ASSM_HEADER 6   E!  8   a   INIT_THERMAL%PRESSURE+SYS_ASSM_HEADER 6   }!  8   a   INIT_THERMAL%VELOCITY+SYS_ASSM_HEADER 7   �!  T   a   SYS_ASSEMBLY%ALLOC+SYS_ASSEMBLY_HEADER 3   	"  N      ALLOC_ASSEMBLY+SYS_ASSEMBLY_HEADER 8   W"  N   a   ALLOC_ASSEMBLY%THIS+SYS_ASSEMBLY_HEADER 7   �"  S   a   SYS_ASSEMBLY%CLEAN+SYS_ASSEMBLY_HEADER 2   �"  k      FREE_ASSEMBLY+SYS_ASSEMBLY_HEADER <   c#  >      FREE_ASSEMBLY%ALLOCATED+SYS_ASSEMBLY_HEADER 7   �#  N   a   FREE_ASSEMBLY%THIS+SYS_ASSEMBLY_HEADER 5   �#  R   a   SYS_ASSEMBLY%SET+SYS_ASSEMBLY_HEADER 1   A$  _      SET_ASSEMBLY+SYS_ASSEMBLY_HEADER 6   �$  N   a   SET_ASSEMBLY%THIS+SYS_ASSEMBLY_HEADER =   �$  N   a   SET_ASSEMBLY%REINPUTDATA+SYS_ASSEMBLY_HEADER 6   <%  S   a   SYS_ASSEMBLY%INIT+SYS_ASSEMBLY_HEADER 2   �%  N      INIT_ASSEMBLY+SYS_ASSEMBLY_HEADER 7   �%  N   a   INIT_ASSEMBLY%THIS+SYS_ASSEMBLY_HEADER 2   +&  N      SET_INPUTDATA+SYS_RE_INPUT_HEADER 7   y&  N   a   SET_INPUTDATA%THIS+SYS_RE_INPUT_HEADER    �&  J       ASSM1 1   '        SYS_RE_INPUT+SYS_RE_INPUT_HEADER 4   $(  @   a   SYS_RE_INPUT%NF+SYS_RE_INPUT_HEADER 4   d(  @   a   SYS_RE_INPUT%NG+SYS_RE_INPUT_HEADER 4   �(  @   a   SYS_RE_INPUT%NS+SYS_RE_INPUT_HEADER 4   �(  @   a   SYS_RE_INPUT%NY+SYS_RE_INPUT_HEADER 6   $)  @   a   SYS_RE_INPUT%NPIN+SYS_RE_INPUT_HEADER 4   d)  @   a   SYS_RE_INPUT%XF+SYS_RE_INPUT_HEADER 4   �)  @   a   SYS_RE_INPUT%XG+SYS_RE_INPUT_HEADER 4   �)  @   a   SYS_RE_INPUT%XS+SYS_RE_INPUT_HEADER 5   $*  @   a   SYS_RE_INPUT%XOS+SYS_RE_INPUT_HEADER 5   d*  @   a   SYS_RE_INPUT%ACF+SYS_RE_INPUT_HEADER 8   �*  @   a   SYS_RE_INPUT%HEIGHT+SYS_RE_INPUT_HEADER 3   �*  @   a   SYS_RE_INPUT%F+SYS_RE_INPUT_HEADER 6   $+  @   a   SYS_RE_INPUT%POUT+SYS_RE_INPUT_HEADER 8   d+  @   a   SYS_RE_INPUT%FLOWIN+SYS_RE_INPUT_HEADER 5   �+  @   a   SYS_RE_INPUT%TIN+SYS_RE_INPUT_HEADER 5   �+  @   a   SYS_RE_INPUT%UIN+SYS_RE_INPUT_HEADER 5   $,  @   a   SYS_RE_INPUT%PIN+SYS_RE_INPUT_HEADER 4   d,  @   a   SYS_RE_INPUT%TI+SYS_RE_INPUT_HEADER 4   �,  @   a   SYS_RE_INPUT%UI+SYS_RE_INPUT_HEADER 4   �,  @   a   SYS_RE_INPUT%PI+SYS_RE_INPUT_HEADER 7   $-  @   a   SYS_RE_INPUT%ALPHA+SYS_RE_INPUT_HEADER 7   d-  @   a   SYS_RE_INPUT%SIGMA+SYS_RE_INPUT_HEADER 5   �-  S   a   SYS_RE_INPUT%SET+SYS_RE_INPUT_HEADER 