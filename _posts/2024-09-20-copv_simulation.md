---
layout: post
title: 'FE Modeling of a COPV using CAM data in Ansys'
date: 2024-09-20
permalink: /posts/2024/09/copv_simulation/
tags:
---
<style>
div {
  text-align: justify;
}

figure {
  text-align: center;
}

figcaption {
  font-style: italic;
  padding: 2px;
  text-align: center;
}

</style>


<br>

<div>

<p>A Composite Overwrapped Pressure Vessel (COPV) has 3 main parts: a metallic boss to facilitate transfer of gas in and out of the vessel, a polymeric liner to prevent diffusion of gas and a composite shell over the liner to provide strength for sustaining the pressure. Filament Winding is the manufacturing process followed to make the composite shell over the liner. To correctly simulate a COPV, it essential to consider the thickness and winding angle variation data from any filament winding software. I have used <a href="https://www.cadfil.com/">CADFIL</a> to get the CAM data and import it in Ansys.</p>

<p>Consider a vessel with radius of 62 mm and length of cylinderical section of 200 mm and both domes are hemispherical. Outer radius of boss opening at one end is 22.5 mm. At the other end, there an adopter of radius of 10 mm attached during manufacturing. Material for Liner is chosen as Polyethylene and Epoxy Carbon UD Wet for composite. The workflow in Ansys Workbench looks like this:</p>

<figure>
  <img  src="/images/copv_simulation/workbench.png" alt="workbench window">
</figure>

<br>

<h3>Making the liner geometry in Ansys</h3>
<div class="two-column-layout">
    <div class="left-column">
        <p>In the Geometry module in Ansys, while drawing a 2D sketch of liner make sure the axis of the cylinder is along the global X axis and the centre of the cylindrical section is at the global origin. This is necessary because the output from Cadfil follows this convention.</p>
    </div>
    <div class="right-column">
            <figure>
            <img  src="/images/copv_simulation/spaceclaim1.png" alt="spaceclaim1" width="150" height="400">
            <img  src="/images/copv_simulation/spaceclaim2.png" alt="spaceclaim2" width="180" height="400">
            <figcaption>Liner Geometry in Spaceclaim</figcaption>
            </figure>
        </div>
</div>


<div class="two-column-layout">
    <div class="left-column">
        <p> Next, make a duplicate of the geometry module and here extract only the outer surface of the liner by selecting the surface and following by Ctrl+c and Ctrl+v. This surface will be used by ACP module to generate the composite layers.</p>
    </div>
    <div class="right-column">
            <figure>
            <img  src="/images/copv_simulation/spaceclaim3.png" alt="spaceclaim3" width="180" height="400">
            <figcaption>Surface for composite layup in Spaceclaim</figcaption>
            </figure>
        </div>
</div>


<div class="two-column-layout">
    <div class="left-column">
        <p> Connect the liner geometry to a 'Mechanical Model' module. Here, assign the material of liner as polyethlyene and mesh the body.</p>
    </div>
    <div class="right-column">
            <figure>
            <img  src="/images/copv_simulation/mechmodel.png" alt="mechmodel">
            <figcaption>Liner Mesh in Mechanical Model</figcaption>
            </figure>
        </div>
</div>


<div class="two-column-layout">
    <div class="left-column">
        <p> Connect the surface geometry to a 'ACP (Pre)' module. Here in the model section, assign the material as Epoxy Carbon and give a random thickness to the surface. After meshing the surface, select the cylindrical surface and add it as a 'Named Selection'. This selection will be used to add hoop layers only in the cylindrical region. Finally, select the circular edges of both the end openings and add them to another Named Selection. These edges will be used as Extrusion Guides for composite layers.</p>
    </div>
    <div class="right-column">
            <figure>
            <img  src="/images/copv_simulation/surfmodel.png" alt="mechmodel">
            <figcaption>Surface with named selection in ACP Pre</figcaption>
            </figure>
        </div>
</div>


<h3>Generating the layup sequence in Cadfil</h3>


<div class="two-column-layout">
    <div class="left-column">
        <p> Before moving forward in Ansys, we need to generate a 'P_lam' file from Cadfil. In the Quickcad menu in Cadfil, select the 'Vessel with endcaps' option. Here, fill the values as desired. In the endcap option, we can choose elliptical(enter height of dome in R2) or torispherical shape (enter knuckle radius in R2 and crown radius in R3). After clicking on 'Calculate', a parameter file (.par) and mandrel file (.mnd) will be created.</p>
    </div>
    <div class="right-column">
            <figure>
            <img  src="/images/copv_simulation/cadfil1.png" alt="cadfil1" width="350" height="450">
            </figure>
        </div>
</div>

<div class="two-column-layout">
    <div class="left-column">
        <p>To set the material data, uncheck the 'use receipe' option and set the values as desired. </p>
    </div>
    <div class="right-column">
            <figure>
            <img  src="/images/copv_simulation/cadfil4.png" alt="cadfil4" width="350" height="400">
            </figure>
        </div>
</div>

<div class="two-column-layout">
    <div class="left-column">
        <p>Next, you will see a band pattern selection table with multiple pattern options. Choose a pattern whose progression factor is closest to one. Other patterns also can be chosen but while actual winding may not be possible because of fibre slippage. </p>
    </div>
    <div class="right-column">
            <figure>
            <img  src="/images/copv_simulation/cadfil2.png" alt="cadfil2">
            </figure>
        </div>
</div>



<div class="two-column-layout">
    <div class="left-column">
        <p>Once everything is done you will see the selected winding pattern and a thickness file (.th2) will be created. Repeat the same procedure for different winding angles and opening radius combination of multiple layers. </p>
    </div>
    <div class="right-column">
            <figure>
            <img  src="/images/copv_simulation/cadfil5.png" alt="cadfil5">
            </figure>
        </div>
</div>


<div class="two-column-layout">
    <div class="left-column">
        <p>Now to generate P_lam file, select 'Analysis(FEA) Output' from Analysis Options menu. Check the Grid point and User points options. Set the name of output file in the first text box. Next choose the mandrel file(.mnd) generated for the first layer. Finally, you have to make the Thickness file list by selecting the .th2 file generated for each layer. The list should be made according to the desired layup sequence. After clicking on calculate, the 'P_lam' file will be generated.  </p>
    </div>
    <div class="right-column">
            <figure>
            <img  src="/images/copv_simulation/cadfil8.png" alt="cadfil8">
            </figure>
        </div>
</div>


<div class="two-column-layout">
    <div class="left-column">
        <p>Now move back to the setup section of ACP Pre. Firstly, in the Units menu set the units to 'mm' because the data from Cadfil is in millimeters. Next, create fabric by selecting the orthotropic material and set a nominal thickness.</p>
    </div>
    <div class="right-column">
            <figure>
            <img  src="/images/copv_simulation/acp1.png" alt="acp1">
            </figure>
        </div>
</div>


<div class="two-column-layout">
    <div class="left-column">
        <p>Next, create an Orientation Selection Set. Select 'All Elements' in the element sets option. This will ensure layers are laid on the entire surface. Then click on the point text box and follow with a click on the surface of the cylinder. Select the "Global Coordinate System' in the Rosettes. Repeat the same procedure for the other element set which we created as a Named Selection for the cylindrical section. This would ensure layers are laid only on the cylindrical surface i.e only hoop layers.</p>
    </div>
    <div class="right-column">
            <figure>
            <img  src="/images/copv_simulation/acp5.png" alt="acp5">
            </figure>
        </div>
</div>

Now, create a 3D Lookup Table and export it to a .csv file. Next, copy the contents of the 'P_lam' generated by Cadfil and paste it in the Lookup Table.csv file generated by ACP as mentioned below. Here, the cells 'column_names', 'column_types' and 'dimensions' must be changed. If there are 10 layers then there would be 20 thickness and 20 angle columns. The column type for both of them is 'scalar'. The dimension for thickness column is 'length' and for the angle column is 'dimensionless'. Once the necessary changes are made, import the modified file back in Ansys ACP.
<br>
<figure>
  <img  src="/images/copv_simulation/plam1.png" alt="plam1">
</figure>

<figure>
  <img  src="/images/copv_simulation/plam2.png" alt="plam2">
</figure>

<figure>
  <img  src="/images/copv_simulation/plam3.png" alt="plam3">
</figure>
<br>

Now in order to make layers in ACP, create a Modeling Group. Next, create a ply and select the orientation set and fabric. Let the ply angle be 0. In the Lookup Table, thickness1, thickness2 and angle1, angle2 contain data for first layer. In the Draping tab, select the tabular values option and in the correction angle 1, select the appropriate column from the Lookup Table. In the thickness tab, again choose From Table option and select the appropriate column. Follow this procedure for all the layers. Hoop layers can be added manually by creating a ply with Orientation Set for cylindrical region.
<br>
<figure>
  <img  src="/images/copv_simulation/acp2.png" alt="acp2" width="320" height="480">
  <img  src="/images/copv_simulation/acp3.png" alt="acp3" width="320" height="480">
  <img  src="/images/copv_simulation/acp4.png" alt="acp4" width="320" height="480">
</figure>
<br>


<div class="two-column-layout">
    <div class="left-column">
        <p>Now to have solid elements for the composite layers, create a solid model and select the options as mentioned in the figure.</p>
    </div>
    <div class="right-column">
            <figure>
            <img  src="/images/copv_simulation/acp6.png" alt="acp6">
            </figure>
        </div>
</div>


<div class="two-column-layout">
    <div class="left-column">
        <p>To provide a guide for extrusion of solid elements, create an extrusion guide for each of the end opening edges we had created in the named selection. In the type option, select direction and provide a unit normal in +X and -X direction in the textbox.</p>
    </div>
    <div class="right-column">
            <figure>
            <img  src="/images/copv_simulation/acp7.png" alt="acp7">
            </figure>
        </div>
</div>


<div class="two-column-layout">
    <div class="left-column">
        <p>After update, the solid model will be generated. We can see the winding angle variation for each layer by clicking on it and turning on 'Show draped fibre angle' option at the top of display window</p>
    </div>
    <div class="right-column">
            <figure>
            <img  src="/images/copv_simulation/acp8.png" alt="acp8">
            </figure>
        </div>
</div>


<div class="two-column-layout">
    <div class="left-column">
        <p>Winding angle for hoop layers can be seen by turning on 'Show fibre direction' option at the top of display window.</p>
    </div>
    <div class="right-column">
            <figure>
            <img  src="/images/copv_simulation/acp9.png" alt="acp9">
            </figure>
        </div>
</div>
<br>

Finally, close the ACP window and add a Static Structural module to workbench. Connect the model of liner to the model section of static structural. And add the setup section of ACP to the same model section of static structural and choose 'Add solid elements' option. After launching the static structural module, you should see the liner and composite assembled together. This can now be used to perform structural simulation as usual.

<figure>
  <img  src="/images/copv_simulation/mech1.png" alt="mech1">
</figure>

<figure>
  <img  src="/images/copv_simulation/mech2.png" alt="mech2">
</figure>

