	  /4  �   k820309    ?          14.0        ШZ                                                                                                           
       E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\debugassembly\driving_output.f90 DRIVING_OUTPUT          RUN_OUTPUT                                                 
                     @               �                'x           #HYDRAU    #GEOM    #MESH %   #PROPERTY 4   #TH_BOUNDARY B   #INITDATA P   #CONFACTOR_ a   #POW j   #THERMAL n   #ALLOC y   #CLEAN |   #SET �   #INIT �   #GRID �            �                                              #HYDRAULIC                  @                               '            #FRIC    #AFLOW    #WET    #DE    #SET 	   #CAL             �                                           	            �                                          	            �                                          	            �                                          	   1     �   � $                      �      	              #SET_HYDRAULIC 
   #     @     @                            
                    #THIS    #FRIC          
                                            #HYDRAULIC          
                                       	  1     �   � $                      �                    #CAL_HYDRAULIC    #     @     @                                                #THIS    #RC    #PD          
                                            #HYDRAULIC          
                                       	        
                                       	           �                                             #ASSMGEOM                  @                               '            #RFUEL    #GASGAP    #SHELLTHICK    #ASSMSHELLTHICK    #ACROSSFLAT    #HEIGHT    #N_PIN    #SET             �                                           	            �                                          	            �                                          	            �                                          	            �                                          	            �                                          	            �                                             1     �   � $                      �                    #SET_ASSMGEOM    #     @     @                                                #THIS    #RFUEL    #GASGAP    #SHELLTHICK     #ASSMSHELLTHICK !   #ACROSSFLAT "   #HEIGHT #   #N_PIN $         
                                            #ASSMGEOM          
                                       	        
                                       	        
                                        	        
                                  !     	        
                                  "     	        
                                  #     	        
                                  $                �                               %     p   ,      #ASSMMESH &                 @              @           &     'p            #NF '   #NG (   #NS )   #NY *   #R +   #Z ,   #SET -            �                               '                        �                               (                       �                               )                       �                               *                     �                               +                 	        &           &                              �                               ,        @         	        &           &                       1     �   � $                      �      -              #SET_ASSMMESH .   #     @     @                            .                    #THIS /   #NF 0   #NG 1   #NS 2   #NY 3         
                                /     p       #ASSMMESH &         
                                  0             
                                  1             
                                  2             
                                  3                �                               4     �   �      #MATERIAL 5                 @              @           5     '�            #RHO 6   #SHC 7   #CTC 8   #DVS 9   #HTC :   #INIT ;          �                               6                  	        &           &                              �                               7        0         	        &           &                              �                               8        `         	        &           &                              �                               9        �         	        &           &                              �                               :        �         	        &                       1     �   � $                      �      ;              #INIT_MATERIAL <   #     @     @                            <                    #THIS =   #NF >   #NG ?   #NS @   #NY A         
                                =     �       #MATERIAL 5         
                                  >             
                                  ?             
                                  @             
                                  A                �                               B        �     #TH_BOUNDARY C                 @                          C     '            #P D   #U H   #T I   #INIT J            �                               D               #BOUNDARY E                 @                          E     '            #INLET F   #OUTLET G            �                               F            	            �                               G           	            �                               H              #BOUNDARY E            �                               I              #BOUNDARY E   1     �   � $                      �      J              #INIT_TH_BOUNDARY K   #     @     @                            K                    #THIS L   #TIN M   #UIN N   #PIN O         
                                L            #TH_BOUNDARY C         
                                  M     	        
                                  N     	        
                                  O     	           �                               P        �     #ASSMINIT Q                 @                          Q     '            #TI R   #PI S   #UI T   #TIN U   #PIN V   #UIN W   #SET X            �                               R            	            �                               S           	            �                               T           	            �                               U           	            �                               V           	            �                               W           	   1     �   � $                      �      X              #SET_ASSMINIT Y   #     @     @                            Y                    #THIS Z   #TI [   #PI \   #UI ]   #TIN ^   #PIN _   #UIN `         
                                Z            #ASSMINIT Q         
                                  [     	        
                                  \     	        
                                  ]     	        
                                  ^     	        
                                  _     	        
                                  `     	           �                               a        �     #CONFACTOR b                 @                          b     '            #ALPHA c   #SIGMA d   #SET e            �                               c            	            �                               d           	   1     �   � $                      �      e              #SET_CONFACTOR f   #     @     @                            f                    #THIS g   #ALPHA h   #SIGMA i         
                                g            #CONFACTOR b         
                                  h     	        
                                  i     	           �                               j     H   �     #ASSMPOW k                 @              @           k     'H            #POWER l   #FQ_CORE m          �                               l                  	        &                              �                               m        $         	        &                                �                               n     x      	   #THERMAL o                 @              @           o     'x            #TEMPERATURE p   #PRESSURE q   #VELOCITY r   #INIT s          �                               p                  	        &           &                              �                               q        0         	        &                              �                               r        T         	        &                       1     �   � $                      �      s              #INIT_THERMAL t   #     @     @                            t                    #THIS u   #TEMPERATURE v   #PRESSURE w   #VELOCITY x         
                                u     x       #THERMAL o         
                                  v     	        
                                  w     	        
                                  x     	  1     �   � $                      �      y         
     #ALLOC_ASSEMBLY z   #     @     @                            z                    #THIS {         
                                {     x      #SYS_ASSEMBLY    1     �   � $                      �      |              #FREE_ASSEMBLY }   #     @     @                            }                   #FREE_ASSEMBLY%ALLOCATED ~   #THIS                                            ~     ALLOCATED       
                                     x      #SYS_ASSEMBLY    1     �   � $                      �      �              #SET_ASSEMBLY �   #     @     @                            �                    #THIS �   #REINPUTDATA �         
                                �     x      #SYS_ASSEMBLY          
                                  �     X      #SYS_RE_INPUT �   1     �   � $                      �      �              #INIT_ASSEMBLY �   #     @     @                            �                    #THIS �         
                                �     x      #SYS_ASSEMBLY    1     �   � $                      �      �              #CAL_GRID �   #     @     @                             �                    #THIS �         
                                �     x      #SYS_ASSEMBLY    #     @                                   �                                   @                           �     'X            #NF �   #NG �   #NS �   #NY �   #NPIN �   #XF �   #XG �   #XS �   #XOS �   #ACF �   #HEIGHT �   #F �   #POUT �   #FLOWIN �   #TIN �   #UIN �   #PIN �   #TI �   #UI �   #PI �   #ALPHA �   #SIGMA �   #SET �            �                               �                        �                               �                       �                               �                       �                               �                       �                               �                       �                               �           	            �                               �           	            �                               �           	            �                               �         	   	            �                               �     $   
   	            �                               �     (      	            �                               �     ,      	            �                               �     0      	            �                               �     4      	            �                               �     8      	            �                               �     <      	            �                               �     @      	            �                               �     D      	            �                               �     H      	            �                               �     L      	            �                               �     P      	            �                               �     T      	   1     �   � $                      �      �              #SET_INPUTDATA �   #     @     @                            �                    #THIS �                                         �     X       #SYS_RE_INPUT �      �   �      fn#fn $   4     b   uapp(DRIVING_OUTPUT    K  <   J  ASSM_GLOBAL 1   �  �       SYS_ASSEMBLY+SYS_ASSEMBLY_HEADER 8   q  O   a   SYS_ASSEMBLY%HYDRAU+SYS_ASSEMBLY_HEADER *   �  |       HYDRAULIC+SYS_ASSM_HEADER /   <  @   a   HYDRAULIC%FRIC+SYS_ASSM_HEADER 0   |  @   a   HYDRAULIC%AFLOW+SYS_ASSM_HEADER .   �  @   a   HYDRAULIC%WET+SYS_ASSM_HEADER -   �  @   a   HYDRAULIC%DE+SYS_ASSM_HEADER .   <  S   a   HYDRAULIC%SET+SYS_ASSM_HEADER .   �  X      SET_HYDRAULIC+SYS_ASSM_HEADER 3   �  K   a   SET_HYDRAULIC%THIS+SYS_ASSM_HEADER 3   2  8   a   SET_HYDRAULIC%FRIC+SYS_ASSM_HEADER .   j  S   a   HYDRAULIC%CAL+SYS_ASSM_HEADER .   �  ^      CAL_HYDRAULIC+SYS_ASSM_HEADER 3     K   a   CAL_HYDRAULIC%THIS+SYS_ASSM_HEADER 1   f  8   a   CAL_HYDRAULIC%RC+SYS_ASSM_HEADER 1   �  8   a   CAL_HYDRAULIC%PD+SYS_ASSM_HEADER 6   �  N   a   SYS_ASSEMBLY%GEOM+SYS_ASSEMBLY_HEADER )   $  �       ASSMGEOM+SYS_ASSM_HEADER /   �  @   a   ASSMGEOM%RFUEL+SYS_ASSM_HEADER 0     @   a   ASSMGEOM%GASGAP+SYS_ASSM_HEADER 4   S  @   a   ASSMGEOM%SHELLTHICK+SYS_ASSM_HEADER 8   �  @   a   ASSMGEOM%ASSMSHELLTHICK+SYS_ASSM_HEADER 4   �  @   a   ASSMGEOM%ACROSSFLAT+SYS_ASSM_HEADER 0   	  @   a   ASSMGEOM%HEIGHT+SYS_ASSM_HEADER /   S	  @   a   ASSMGEOM%N_PIN+SYS_ASSM_HEADER -   �	  R   a   ASSMGEOM%SET+SYS_ASSM_HEADER -   �	  �      SET_ASSMGEOM+SYS_ASSM_HEADER 2   �
  J   a   SET_ASSMGEOM%THIS+SYS_ASSM_HEADER 3   �
  8   a   SET_ASSMGEOM%RFUEL+SYS_ASSM_HEADER 4     8   a   SET_ASSMGEOM%GASGAP+SYS_ASSM_HEADER 8   O  8   a   SET_ASSMGEOM%SHELLTHICK+SYS_ASSM_HEADER <   �  8   a   SET_ASSMGEOM%ASSMSHELLTHICK+SYS_ASSM_HEADER 8   �  8   a   SET_ASSMGEOM%ACROSSFLAT+SYS_ASSM_HEADER 4   �  8   a   SET_ASSMGEOM%HEIGHT+SYS_ASSM_HEADER 3   /  8   a   SET_ASSMGEOM%N_PIN+SYS_ASSM_HEADER 6   g  N   a   SYS_ASSEMBLY%MESH+SYS_ASSEMBLY_HEADER )   �  {       ASSMMESH+SYS_ASSM_HEADER ,   0  @   a   ASSMMESH%NF+SYS_ASSM_HEADER ,   p  @   a   ASSMMESH%NG+SYS_ASSM_HEADER ,   �  @   a   ASSMMESH%NS+SYS_ASSM_HEADER ,   �  @   a   ASSMMESH%NY+SYS_ASSM_HEADER +   0  |   a   ASSMMESH%R+SYS_ASSM_HEADER +   �  |   a   ASSMMESH%Z+SYS_ASSM_HEADER -   (  R   a   ASSMMESH%SET+SYS_ASSM_HEADER -   z  n      SET_ASSMMESH+SYS_ASSM_HEADER 2   �  J   a   SET_ASSMMESH%THIS+SYS_ASSM_HEADER 0   2  8   a   SET_ASSMMESH%NF+SYS_ASSM_HEADER 0   j  8   a   SET_ASSMMESH%NG+SYS_ASSM_HEADER 0   �  8   a   SET_ASSMMESH%NS+SYS_ASSM_HEADER 0   �  8   a   SET_ASSMMESH%NY+SYS_ASSM_HEADER :     N   a   SYS_ASSEMBLY%PROPERTY+SYS_ASSEMBLY_HEADER )   `  {       MATERIAL+SYS_ASSM_HEADER -   �  |   a   MATERIAL%RHO+SYS_ASSM_HEADER -   W  |   a   MATERIAL%SHC+SYS_ASSM_HEADER -   �  |   a   MATERIAL%CTC+SYS_ASSM_HEADER -   O  |   a   MATERIAL%DVS+SYS_ASSM_HEADER -   �  l   a   MATERIAL%HTC+SYS_ASSM_HEADER .   7  S   a   MATERIAL%INIT+SYS_ASSM_HEADER .   �  n      INIT_MATERIAL+SYS_ASSM_HEADER 3   �  J   a   INIT_MATERIAL%THIS+SYS_ASSM_HEADER 1   B  8   a   INIT_MATERIAL%NF+SYS_ASSM_HEADER 1   z  8   a   INIT_MATERIAL%NG+SYS_ASSM_HEADER 1   �  8   a   INIT_MATERIAL%NS+SYS_ASSM_HEADER 1   �  8   a   INIT_MATERIAL%NY+SYS_ASSM_HEADER =   "  Q   a   SYS_ASSEMBLY%TH_BOUNDARY+SYS_ASSEMBLY_HEADER ,   s  c       TH_BOUNDARY+SYS_ASSM_HEADER .   �  N   a   TH_BOUNDARY%P+SYS_ASSM_HEADER )   $  [       BOUNDARY+SYS_ASSM_HEADER /     @   a   BOUNDARY%INLET+SYS_ASSM_HEADER 0   �  @   a   BOUNDARY%OUTLET+SYS_ASSM_HEADER .   �  N   a   TH_BOUNDARY%U+SYS_ASSM_HEADER .   M  N   a   TH_BOUNDARY%T+SYS_ASSM_HEADER 1   �  V   a   TH_BOUNDARY%INIT+SYS_ASSM_HEADER 1   �  i      INIT_TH_BOUNDARY+SYS_ASSM_HEADER 6   Z  M   a   INIT_TH_BOUNDARY%THIS+SYS_ASSM_HEADER 5   �  8   a   INIT_TH_BOUNDARY%TIN+SYS_ASSM_HEADER 5   �  8   a   INIT_TH_BOUNDARY%UIN+SYS_ASSM_HEADER 5     8   a   INIT_TH_BOUNDARY%PIN+SYS_ASSM_HEADER :   O  N   a   SYS_ASSEMBLY%INITDATA+SYS_ASSEMBLY_HEADER )   �  �       ASSMINIT+SYS_ASSM_HEADER ,     @   a   ASSMINIT%TI+SYS_ASSM_HEADER ,   ]  @   a   ASSMINIT%PI+SYS_ASSM_HEADER ,   �  @   a   ASSMINIT%UI+SYS_ASSM_HEADER -   �  @   a   ASSMINIT%TIN+SYS_ASSM_HEADER -     @   a   ASSMINIT%PIN+SYS_ASSM_HEADER -   ]  @   a   ASSMINIT%UIN+SYS_ASSM_HEADER -   �  R   a   ASSMINIT%SET+SYS_ASSM_HEADER -   �  �      SET_ASSMINIT+SYS_ASSM_HEADER 2   p  J   a   SET_ASSMINIT%THIS+SYS_ASSM_HEADER 0   �  8   a   SET_ASSMINIT%TI+SYS_ASSM_HEADER 0   �  8   a   SET_ASSMINIT%PI+SYS_ASSM_HEADER 0   *  8   a   SET_ASSMINIT%UI+SYS_ASSM_HEADER 1   b  8   a   SET_ASSMINIT%TIN+SYS_ASSM_HEADER 1   �  8   a   SET_ASSMINIT%PIN+SYS_ASSM_HEADER 1   �  8   a   SET_ASSMINIT%UIN+SYS_ASSM_HEADER <   
  O   a   SYS_ASSEMBLY%CONFACTOR_+SYS_ASSEMBLY_HEADER *   Y  c       CONFACTOR+SYS_ASSM_HEADER 0   �  @   a   CONFACTOR%ALPHA+SYS_ASSM_HEADER 0   �  @   a   CONFACTOR%SIGMA+SYS_ASSM_HEADER .   <   S   a   CONFACTOR%SET+SYS_ASSM_HEADER .   �   d      SET_CONFACTOR+SYS_ASSM_HEADER 3   �   K   a   SET_CONFACTOR%THIS+SYS_ASSM_HEADER 4   >!  8   a   SET_CONFACTOR%ALPHA+SYS_ASSM_HEADER 4   v!  8   a   SET_CONFACTOR%SIGMA+SYS_ASSM_HEADER 5   �!  M   a   SYS_ASSEMBLY%POW+SYS_ASSEMBLY_HEADER (   �!  \       ASSMPOW+SYS_ASSM_HEADER .   W"  l   a   ASSMPOW%POWER+SYS_ASSM_HEADER 0   �"  l   a   ASSMPOW%FQ_CORE+SYS_ASSM_HEADER 9   /#  M   a   SYS_ASSEMBLY%THERMAL+SYS_ASSEMBLY_HEADER (   |#  {       THERMAL+SYS_ASSM_HEADER 4   �#  |   a   THERMAL%TEMPERATURE+SYS_ASSM_HEADER 1   s$  l   a   THERMAL%PRESSURE+SYS_ASSM_HEADER 1   �$  l   a   THERMAL%VELOCITY+SYS_ASSM_HEADER -   K%  R   a   THERMAL%INIT+SYS_ASSM_HEADER -   �%  {      INIT_THERMAL+SYS_ASSM_HEADER 2   &  I   a   INIT_THERMAL%THIS+SYS_ASSM_HEADER 9   a&  8   a   INIT_THERMAL%TEMPERATURE+SYS_ASSM_HEADER 6   �&  8   a   INIT_THERMAL%PRESSURE+SYS_ASSM_HEADER 6   �&  8   a   INIT_THERMAL%VELOCITY+SYS_ASSM_HEADER 7   	'  T   a   SYS_ASSEMBLY%ALLOC+SYS_ASSEMBLY_HEADER 3   ]'  N      ALLOC_ASSEMBLY+SYS_ASSEMBLY_HEADER 8   �'  N   a   ALLOC_ASSEMBLY%THIS+SYS_ASSEMBLY_HEADER 7   �'  S   a   SYS_ASSEMBLY%CLEAN+SYS_ASSEMBLY_HEADER 2   L(  k      FREE_ASSEMBLY+SYS_ASSEMBLY_HEADER <   �(  >      FREE_ASSEMBLY%ALLOCATED+SYS_ASSEMBLY_HEADER 7   �(  N   a   FREE_ASSEMBLY%THIS+SYS_ASSEMBLY_HEADER 5   C)  R   a   SYS_ASSEMBLY%SET+SYS_ASSEMBLY_HEADER 1   �)  _      SET_ASSEMBLY+SYS_ASSEMBLY_HEADER 6   �)  N   a   SET_ASSEMBLY%THIS+SYS_ASSEMBLY_HEADER =   B*  N   a   SET_ASSEMBLY%REINPUTDATA+SYS_ASSEMBLY_HEADER 6   �*  S   a   SYS_ASSEMBLY%INIT+SYS_ASSEMBLY_HEADER 2   �*  N      INIT_ASSEMBLY+SYS_ASSEMBLY_HEADER 7   1+  N   a   INIT_ASSEMBLY%THIS+SYS_ASSEMBLY_HEADER 6   +  N   a   SYS_ASSEMBLY%GRID+SYS_ASSEMBLY_HEADER -   �+  N      CAL_GRID+SYS_ASSEMBLY_HEADER 2   ,  N   a   CAL_GRID%THIS+SYS_ASSEMBLY_HEADER    i,  D       RUN_OUTPUT 1   �,        SYS_RE_INPUT+SYS_RE_INPUT_HEADER 4   �-  @   a   SYS_RE_INPUT%NF+SYS_RE_INPUT_HEADER 4    .  @   a   SYS_RE_INPUT%NG+SYS_RE_INPUT_HEADER 4   @.  @   a   SYS_RE_INPUT%NS+SYS_RE_INPUT_HEADER 4   �.  @   a   SYS_RE_INPUT%NY+SYS_RE_INPUT_HEADER 6   �.  @   a   SYS_RE_INPUT%NPIN+SYS_RE_INPUT_HEADER 4    /  @   a   SYS_RE_INPUT%XF+SYS_RE_INPUT_HEADER 4   @/  @   a   SYS_RE_INPUT%XG+SYS_RE_INPUT_HEADER 4   �/  @   a   SYS_RE_INPUT%XS+SYS_RE_INPUT_HEADER 5   �/  @   a   SYS_RE_INPUT%XOS+SYS_RE_INPUT_HEADER 5    0  @   a   SYS_RE_INPUT%ACF+SYS_RE_INPUT_HEADER 8   @0  @   a   SYS_RE_INPUT%HEIGHT+SYS_RE_INPUT_HEADER 3   �0  @   a   SYS_RE_INPUT%F+SYS_RE_INPUT_HEADER 6   �0  @   a   SYS_RE_INPUT%POUT+SYS_RE_INPUT_HEADER 8    1  @   a   SYS_RE_INPUT%FLOWIN+SYS_RE_INPUT_HEADER 5   @1  @   a   SYS_RE_INPUT%TIN+SYS_RE_INPUT_HEADER 5   �1  @   a   SYS_RE_INPUT%UIN+SYS_RE_INPUT_HEADER 5   �1  @   a   SYS_RE_INPUT%PIN+SYS_RE_INPUT_HEADER 4    2  @   a   SYS_RE_INPUT%TI+SYS_RE_INPUT_HEADER 4   @2  @   a   SYS_RE_INPUT%UI+SYS_RE_INPUT_HEADER 4   �2  @   a   SYS_RE_INPUT%PI+SYS_RE_INPUT_HEADER 7   �2  @   a   SYS_RE_INPUT%ALPHA+SYS_RE_INPUT_HEADER 7    3  @   a   SYS_RE_INPUT%SIGMA+SYS_RE_INPUT_HEADER 5   @3  S   a   SYS_RE_INPUT%SET+SYS_RE_INPUT_HEADER 2   �3  N      SET_INPUTDATA+SYS_RE_INPUT_HEADER 7   �3  N   a   SET_INPUTDATA%THIS+SYS_RE_INPUT_HEADER 