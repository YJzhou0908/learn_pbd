## Mass-Spring 

Result:

<div style="text-align: center;">
<img src="./gif/xpbd_spring.gif" style="width: 30%; height: auto;" />
</div>


## FEM(Based on Neo-Hookean Hyperelasticity Model)

Hydrostatic Constraint:

$$
C_H(\mathbf{F}) = \det(\mathbf{F}) - 1
$$

Deviatoric Constraint:

$$
C_D(\mathbf{F}) = \sqrt{\text{tr}(\mathbf{F}^T \mathbf{F})}
$$

Result:

*If the gradient of hydrostatic constraint is wrong:*
<div style="text-align: center;">
<img src="./gif/xpbd_ball1.gif" style="width: 50%; height: auto;" />
</div>

*Right simulation:*
<div style="text-align: center;">
<img src="./gif/xpbd_right.gif" style="width: 50%; height: auto;" />
</div>

