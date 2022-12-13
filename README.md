# dorsal ventral analysis
 This is the library used to align dorsal and ventrla view of a wing. It also generates the second-level product, whcih contains grided multispectral reflectance data, wing shapes, and measurements of body size and antenna.
 
Please find the article describing the application of this script below (DOI: 10.1038/s42003-022-04282-z):
(https://www.nature.com/articles/s42003-022-04282-z)

Exemplar datasets can be found in this repository:
(https://doi.org/10.5061/dryad.37pvmcvp5)

Please find the link for the corresponding protocol:
(https://www.protocols.io/private/F3292DF1FE0211EB878B0A58A9FEAC02)

Please cite this work as:
Chan, W-P., Rabideau Childers, R., Ashe, S. *et al.* A high-throughput multispectral imaging system for museum specimens. *Commun Biol* **5**, 1318 (2022).

 
## Output data sturecutre (also included in the corresponding scripts):
      %{1}=All parameters and grids
            %{1}{1}=dorsal
            %{1}{2}=ventral
                % {1}{}{1}={seg4PtsLF,wingGridsLF }; %4 key points & grids
                % {1}{}{2}={seg4PtsRF,wingGridsRF}; %4 key points & grids
                % {1}{}{3}={outSeg4PtsLH ,wingGridsLH}; %4 key points & grids
                % {1}{}{4}={outSeg4PtsRH ,wingGridsRH}; %4 key points & grids
                % {1}{}{5}={refineAreaLH,regPtLH,reconstructRegPtLH}; %mask without ornament, 3 key points & key points of the mask without ornament
                % {1}{}{6}={refineAreaRH,regPtRH,reconstructRegPtRH}; %mask without ornament, 3 key points & key points of the mask without ornament
                % {1}{}{7}=original mask of LH
                % {1}{}{8}=original mask of RH
     %{2}=All reflectance grids summary; if it's not applicable, use -9999 as a placeholder
            %{2}{1}=dorsal
            %{2}{2}=ventral
                %{2}{}{1}=740
                %{2}{}{2}=940
                %{2}{}{3}=UV
                %{2}{}{4}=UVF
                %{2}{}{5}=F
                %{2}{}{6}=white
                %{2}{}{7}=whitePo1
                %{2}{}{8}=whitePo2
                %{2}{}{9}=FinRGB
                %{2}{}{10}=PolDiff
                    %{2}{}{}{1}=left forewing
                    %{2}{}{}{2}=right forewing
                    %{2}{}{}{3}=left hindwing
                    %{2}{}{}{4}=right hindwing
                        %{2}{}{}{}{1}=mean reflectance in a grid (gray scale: single layer; RGB: 3 layers)
                        %{2}{}{}{}{2}=standard error of reflectance in a grid (gray scale: single layer; RGB: 3 layers)
      %{3}=scale length of dorsal and ventral images
            %{3}(1)=dorsal
            %{3}(2)=ventral
      %{4}=original file information
            %{4}(1)=dorsal side, morph-seg file name
            %{4}(2)=ventral side, morph-seg file name
            %{4}(3)=dorsal side, AllBandsMask file name
            %{4}(4)=ventral side, AllBandsMask file name
        %{5}=tails information
            %{5}{1}=dorsal_LH tails (-9999 represents no tail)
            %{5}{2}=ventral_LH tails (-9999 represents no tail)
            %{5}{3}=dorsal_RH tails (-9999 represents no tail)
            %{5}{4}=ventral_RH tails (-9999 represents no tail)
                %{5}{}{N}= No. N tail part
                    %{5}{}{}{1}= tail mask
                    %{5}{}{}{2}= raw coordinations of tail base (the part connect to hindwing; points of two ends and the mid-point are provided [mid-point is at second row])
                    %{5}{}{}{3}= summary grid coordinations of tail base (the part connect to hindwing; points of two ends and the mid-point are provided [mid-point is at second row])
                    %{5}{}{}{4}= [tail-base boundary Length, tail length, tail area, tail width, tailCurvature] (in cm; curvature is unit less)
        %{6}= information related to body size and antennae
            %{6}{1}=body length and width (in cm) (dorsal and ventral sides are in two rows)
            %{6}{2}= Antenna length, width, bulb width, degree of curvature (all in mm); LEFT antenna; first row is dorsal side, second is ventral side.
            %{6}{3}= Antenna length, width, bulb width, degree of curvature (all in mm); RIGHT antenna; first row is dorsal side, second is ventral side.
        %{7}= All reflectance on antennae; if it's not applicable, use -9999 as a placeholder
            %{7}{1}= Mean reflectance rescale to 100 grids from the tip of antenna to its base
            %{7}{2}= Mean reflectance every 0.1 mm from the tip of antenna to its base
                %{7}{1|2}{1}= Left antenna
                %{7}{1|2}{2}= Right antenna
                    %{7}{1|2}{}(: , : , 1)= dorsal side
                    %{7}{1|2}{}(: , : , 2)= ventral side
                    %{7}{1|2}{}(: , 1 , :)= distance to the tip of antenna (mm)
                    %{7}{1|2}{}(: , 2:21, :)= spectral reflectance (740, 940, UV, UVF, F, white (R, G, B), whitePo1 (R, G, B), whitePo2 (R, G, B), FinRGB (R, G, B), PolDiff (R, G, B))
            %{7}{3}= Raw reflectance value (raw pixel-based value)
                %{7}{3}{1}= Left antenna
                %{7}{3}{2}= Right antenna
                    %{7}{3}{}{1}= reflectance matrix
                        %{7}{3}{}{1}(:,1)= Percentage from the tip
                        %{7}{3}{}{1}(:,2:21)= spectral reflectance (740, 940, UV, UVF, F, white (R, G, B), whitePo1 (R, G, B), whitePo2 (R, G, B), FinRGB (R, G, B), PolDiff (R, G, B))
                        %{7}{3}{}{1}(:,22)= pixel distance from the tip
                        %{7}{3}{}{1}(:,23)= original coordination X
                        %{7}{3}{}{1}(:,24)= original coordination Y
                    %{7}{3}{}{2}= scale (= 1 cm)

