from math import *
from manim import *

#manim -pqk scene.py Part

# Colors
WHITE = '#FFFFFF'
BLACK = '#000000'
RED = '#C02020'
GREEN = '#70FF70'
BLUE = '#8080FF'
PURE_RED = '#FF0000'
PURE_GREEN = '#00FF00'
PURE_BLUE = '#0000FF'

# Unit & ORIGIN vectors
# ORIGIN ; UP ; DOWN ; LEFT ; RIGHT ; IN ; OUT

class Part2c_4D(ThreeDScene):
    def Lerp(self , starting_point , ending_point , t , flag=0):
        return [int(a*(1-t)+b*t) for a, b in zip(starting_point, ending_point)] if flag else [a*(1-t)+b*t for a, b in zip(starting_point, ending_point)]
    def InverseLerp(self, starting_point , ending_point, x):
        return (x-starting_point)/(ending_point-starting_point)
    
    def construct(self):
        self.camera.background_color = '#101020'
        self.camera.light_source.move_to(OUT*20)
        self.set_camera_orientation(phi=2*PI/5, theta=PI/4)
        self.set_camera_orientation(zoom=1.5)
        axes = ThreeDAxes(
            x_range=(-3, 3, 1), y_range=(-3, 3, 1), z_range=(-2, 2, 1),
            x_length=6, y_length=6, z_length=4
        ).set_opacity(0.25)
        x_axis_label_rotation_angle , x_axis_label_rotation_axis = angle_axis_from_quaternion(quaternion_mult(quaternion_from_angle_axis(angle=PI, axis=OUT, axis_normalized=True) , quaternion_from_angle_axis(angle=PI/2, axis=RIGHT, axis_normalized=True)))
        y_axis_label_rotation_angle , y_axis_label_rotation_axis = angle_axis_from_quaternion(quaternion_mult(quaternion_from_angle_axis(angle=PI/2, axis=OUT, axis_normalized=True) , quaternion_from_angle_axis(angle=PI/2, axis=RIGHT, axis_normalized=True)))
        x_axis_label = axes.get_x_axis_label(MathTex("x" , color = RED).scale(1) , edge=RIGHT , buff=1).move_to([3.25,0,0]).rotate(x_axis_label_rotation_angle , axis=x_axis_label_rotation_axis)
        y_axis_label = axes.get_y_axis_label(MathTex("y" , color = GREEN).scale(1) , edge=UP , buff=1 , rotation=0).move_to([0,3.25,0]).rotate(y_axis_label_rotation_angle , axis=y_axis_label_rotation_axis)
        z_axis_label = axes.get_z_axis_label(MathTex("z" , color = BLUE).scale(1) , edge=OUT , buff=1 , rotation=PI/2).move_to([0,0,2.25]).rotate(PI*3/4 , axis=OUT)
        axes_labels = Group(x_axis_label , y_axis_label , z_axis_label)

        #define t
        t_label = Variable(var=0.00, label=Text('t'), num_decimal_places=2).scale(1.5).to_corner(UP + RIGHT)
        t_label.label.set_color(GREEN)
        t_label.value.set_color(BLUE)
        # define f1 and f2
        f1 = np.array([0.29,0.00,0.03,0.68])
        f2 = np.array([0.10,0.09,0.12,0.69])
        
        # create a
        a = always_redraw(lambda : Line3D(ORIGIN, self.Lerp(f1,f2,t_label.value.get_value())[:3], color=ManimColor.from_rgb(self.Lerp([0,0,255],[255,0,0],self.InverseLerp(0.755,0.831,self.Lerp(f1[3:],f2[3:],t_label.value.get_value())[0])),1), thickness=0.01))

        # labels
        self.add_fixed_in_frame_mobjects(t_label.label)
        self.add_fixed_in_frame_mobjects(t_label.value)

        # fade in scene
        self.play(FadeIn(axes) , FadeIn(axes_labels) , FadeIn(t_label.label) , FadeIn(t_label.value) , run_time=2)
        self.add(t_label.label, t_label.value)
        self.wait(0.5)
        self.play(GrowFromPoint(a,ORIGIN))
        self.wait(1)
































class Part2c_Lerp(Scene):
    def Lerp(self , starting_point , ending_point , t , flag=0):
        return [int(a*(1-t)+b*t) for a, b in zip(starting_point, ending_point)] if flag else [a*(1-t)+b*t for a, b in zip(starting_point, ending_point)]
    
    def CircleCoordinate(self , starting_point , ending_point , t):
        x = MathTex(r"[\:" , str(round(starting_point[0]*(1-t)+ending_point[0]*t , 1)) , r"\:,\;" , str(-2.0) , r"\:]").set_color(color=RED)
        x[1].set_color(color=WHITE)
        return x
    
    def CircleColor(self , starting_point , ending_point , t):
        color = self.Lerp(starting_point , ending_point , t , 1)
        a = len(str(color[0]))
        b = len(str(color[1]))
        c = len(str(color[2]))
        x = MathTex(r"\begin{bmatrix}"+str(color[0])+r"\\"+str(color[1])+r"\\"+str(color[2])+r"\end{bmatrix}" , font_size=96)
        for i in range(2,2+a):
            x[0][i].set_color(color='PURE_RED')
        for j in range(2+a,2+a+b):
            x[0][j].set_color(color='PURE_GREEN')
        for k in range(2+a+b,2+a+b+c):
            x[0][k].set_color(color='PURE_BLUE')
        return x

    def construct(self):
        self.camera.background_color = '#101028'
        #Let's revisit the 3 examples from the beginning.
        t_label = Variable(var=0.00, label=Text('t'), num_decimal_places=2).scale(1.5).move_to([0,3,0])
        t_label.label.set_color(GREEN)
        t_label.value.set_color(BLUE)
        starting_point = [-3,-2,0]
        ending_point = [3,-2,0]
        circle = always_redraw(lambda : Circle(arc_center=self.Lerp(starting_point , ending_point , t_label.value.get_value()) , radius=1 , fill_color=BLACK , fill_opacity=1 , stroke_color=RED, stroke_width=20))
        circle_coordinate = always_redraw(lambda : self.CircleCoordinate(starting_point , ending_point , t_label.value.get_value()).next_to(circle , UP))

        lerp_equation = MathTex(r"(1-" , r"t" , r") \cdot" , r"s_{tart}" , r"+" , r"t" , r"\cdot" , r"e_{nd}" , font_size=72 , color=WHITE).move_to([0,1.5,0])
        lerp_equation[1].set_fill(color=GREEN)
        lerp_equation[5].set_fill(color=GREEN)
        lerp_equation[3].set_fill(color=RED)
        lerp_equation[7].set_fill(color=RED)

        self.wait(0.5)
        self.play(Create(circle) , FadeIn(t_label) , FadeIn(lerp_equation) , run_time=2)
        self.play(Wait(run_time=1))

#When animating movement, a is set as the initial coordinate, and b is the final coordinate.
        lerp_equation.generate_target()
        lerp_equation.target = MathTex(r"(1-" , r"t" , r") \cdot" , r"\begin{bmatrix} -3 \\ -2 \end{bmatrix}" , r"+" , r"t" , r"\cdot" , r"e_{nd}" , font_size=72 , color=WHITE).move_to([0,1.5,0])
        lerp_equation.target[1].set_fill(color=GREEN)
        lerp_equation.target[5].set_fill(color=GREEN)
        lerp_equation.target[3].set_fill(color=RED)
        lerp_equation.target[7].set_fill(color=RED)
        self.play(MoveToTarget(lerp_equation) , run_time=0.5)

        self.play(Wait(run_time=1))
        lerp_equation.generate_target()
        lerp_equation.target = MathTex(r"(1-" , r"t" , r") \cdot" , r"\begin{bmatrix} -3 \\ -2 \end{bmatrix}" , r"+" , r"t" , r"\cdot" , r"\begin{bmatrix} 3 \\ -2 \end{bmatrix}" , font_size=72 , color=WHITE).move_to([0,1.5,0])
        lerp_equation.target[1].set_fill(color=GREEN)
        lerp_equation.target[5].set_fill(color=GREEN)
        lerp_equation.target[3].set_fill(color=RED)
        lerp_equation.target[7].set_fill(color=RED)
        self.play(MoveToTarget(lerp_equation) , run_time=0.5)
        self.play(Wait(run_time=1))

#Since a and b are in the form of coordinates, or vectors, the output of this formula will also be a coordinate.
        self.play(FadeIn(circle_coordinate) , run_time=1)
        self.play(Wait(run_time=1))

#As t increases from 0 to 1, the object moves from its starting position to its final position
        self.play(t_label.tracker.animate.set_value(1) , run_time=7)
        self.play(t_label.tracker.animate.set_value(1) , run_time=2)

#Lerp is also used to blend colors
        self.clear()
        starting_point=[160 , 0 , 0]
        ending_point=[15 , 255 , 15]
        t_label.tracker.set_value(0)
        self.add(t_label)
        circle = always_redraw(lambda : Square(side_length=3 , fill_color=ManimColor.from_rgb(self.Lerp(starting_point , ending_point , t_label.value.get_value() , 1)), fill_opacity=1 , stroke_color=BLACK , stroke_width=10).move_to([-2,-2,0]))
        self.add(circle)
        lerp_equation = MathTex(r"(1-" , r"t" , r") \cdot" , r"s_{tart}" , r"+" , r"t" , r"\cdot" , r"e_{nd}" , font_size=54 , color=WHITE).move_to([0,1.25,0])
        lerp_equation[1].set_fill(color=GREEN)
        lerp_equation[5].set_fill(color=GREEN)
        lerp_equation[3].set_fill(color=RED)
        lerp_equation[7].set_fill(color=RED)
        self.add(lerp_equation)
        self.play(t_label.tracker.animate.set_value(0) , run_time=1)

#In this case, a is the RGB value of the initial color expressed as a vector, and b is the RGB value of the final color
        self.play(Wait(run_time=1))
        lerp_equation.generate_target()
        lerp_equation.target = MathTex(r"(1-" , r"t" , r") \cdot" , r"\begin{bmatrix} 160 \\ 0 \\ 0 \end{bmatrix}" , r"+" , r"t" , r"\cdot" , r"e_{nd}" , font_size=54 , color=WHITE).move_to([0,1.25,0])
        lerp_equation.target[1].set_fill(color=GREEN)
        lerp_equation.target[5].set_fill(color=GREEN)
        lerp_equation.target[3].set_fill(color=RED)
        lerp_equation.target[3][2].set_fill(color='PURE_RED')
        lerp_equation.target[3][3].set_fill(color='PURE_RED')
        lerp_equation.target[3][4].set_fill(color='PURE_RED')
        lerp_equation.target[3][5].set_fill(color='PURE_GREEN')
        lerp_equation.target[3][6].set_fill(color='#0040FF')
        lerp_equation.target[7].set_fill(color=RED)
        self.play(MoveToTarget(lerp_equation) , run_time=0.5)

        self.play(Wait(run_time=1))
        lerp_equation.generate_target()
        lerp_equation.target = MathTex(r"(1-" , r"t" , r") \cdot" , r"\begin{bmatrix} 160 \\ 0 \\ 0 \end{bmatrix}" , r"+" , r"t" , r"\cdot" , r"\begin{bmatrix} 15 \\ 225 \\ 15 \end{bmatrix}" , font_size=54 , color=WHITE).move_to([0,1.25,0])
        lerp_equation.target[1].set_fill(color=GREEN)
        lerp_equation.target[5].set_fill(color=GREEN)
        lerp_equation.target[3].set_fill(color=RED)
        lerp_equation.target[3][2].set_fill(color='PURE_RED')
        lerp_equation.target[3][3].set_fill(color='PURE_RED')
        lerp_equation.target[3][4].set_fill(color='PURE_RED')
        lerp_equation.target[3][5].set_fill(color='PURE_GREEN')
        lerp_equation.target[3][6].set_fill(color='#0050FF')
        lerp_equation.target[7].set_fill(color=RED)
        lerp_equation.target[7][2].set_fill(color='PURE_RED')
        lerp_equation.target[7][3].set_fill(color='PURE_RED')
        lerp_equation.target[7][4].set_fill(color='PURE_GREEN')
        lerp_equation.target[7][5].set_fill(color='PURE_GREEN')
        lerp_equation.target[7][6].set_fill(color='PURE_GREEN')
        lerp_equation.target[7][7].set_fill(color='#0050FF')
        lerp_equation.target[7][8].set_fill(color='#0050FF')
        self.play(MoveToTarget(lerp_equation) , run_time=0.5)
        self.play(Wait(run_time=1))

#This way, the output of this equation will also be an RGB colorcode.
        circle_color = always_redraw(lambda : self.CircleColor(starting_point , ending_point , t_label.value.get_value()).next_to(circle , RIGHT*3))
        self.play(FadeIn(circle_color) , run_time=2)
        self.play(Wait(run_time=2))

#As t increases from 0 to 1, the object's color changes from its starting color to its final color
        self.play(t_label.tracker.animate.set_value(1) , run_time=6)
        self.play(Wait(run_time=2))

#The last example in the beginning was transforming between two graphs
        self.clear()
        coordinate_plane = NumberPlane(background_line_style={
                "stroke_color": BLUE_D,
                "stroke_width": 2,
                "stroke_opacity": 0.5
            }).add_coordinates().set_opacity(0.5)
        self.add(coordinate_plane)
        t_label.move_to([0 , 3 , 0]).set_z_index(3)
        t_label.tracker.set_value(0)
        t_label.value.set_value(0)
        self.add(t_label)
        f_x = lambda x : (x**2)*sin(x)
        g_x = lambda x : cos(x)
        h_x = lambda x : (1-t_label.value.get_value())*f_x(x)  + (t_label.value.get_value())*g_x(x)
        graph_f = FunctionGraph(f_x , color=BLUE).set_stroke(opacity=0.3).set_z_index(1)
        graph_g = FunctionGraph(g_x , color=BLUE).set_stroke(opacity=0.3).set_z_index(1)
        graph_h = always_redraw(lambda : FunctionGraph(h_x , color='#F0F050').set_z_index(2))
        coordinate_plane.add(graph_f)
        coordinate_plane.add(graph_g)
        coordinate_plane.add(graph_h)
        lerp_equation = MathTex(r"(1-" , r"t" , r") \cdot" , r"s_{tart}" , r"+" , r"t" , r"\cdot" , r"e_{nd}" , font_size=72 , color=WHITE).move_to([0 , 2 , 0])
        lerp_equation[1].set_fill(color=GREEN)
        lerp_equation[5].set_fill(color=GREEN)
        lerp_equation[3].set_fill(color='#FF5050')
        lerp_equation[7].set_fill(color='#FF5050')
        lerp_equation.set_z_index(4)
        self.add(lerp_equation)
        self.play(Wait(run_time=2) , t_label.tracker.animate.set_value(0))

#Here, a is the initial function, and b is the final function.
        lerp_equation.generate_target()
        lerp_equation.target = MathTex(r"(1-" , r"t" , r") \cdot" , r"f(x)" , r"+" , r"t" , r"\cdot" , r"e_{nd}" , font_size=72 , color=WHITE).move_to([0,2,0])
        lerp_equation.target[1].set_fill(color=GREEN)
        lerp_equation.target[5].set_fill(color=GREEN)
        lerp_equation.target[3].set_fill(color='#FF5050')
        lerp_equation.target[7].set_fill(color='#FF5050')
        lerp_equation.set_z_index(4)
        self.play(MoveToTarget(lerp_equation) , run_time=0.5)

        self.play(Wait(run_time=1))
        lerp_equation.generate_target()
        lerp_equation.target = MathTex(r"(1-" , r"t" , r") \cdot" , r"f(x)" , r"+" , r"t" , r"\cdot" , r"g(x)" , font_size=72 , color=WHITE).move_to([0,2,0])
        lerp_equation.target[1].set_fill(color=GREEN)
        lerp_equation.target[5].set_fill(color=GREEN)
        lerp_equation.target[3].set_fill(color='#FF5050')
        lerp_equation.target[7].set_fill(color='#FF5050')
        lerp_equation.set_z_index(4)
        self.play(MoveToTarget(lerp_equation) , run_time=0.5)
        self.play(Wait(run_time=1))

#As t increases, the output gradually transforms from f(x) to g(x)
        self.play(t_label.tracker.animate.set_value(1) , run_time=6 , func_rate=smooth)
        self.play(t_label.tracker.animate.set_value(1) , run_time=2)
        cover = Rectangle(stroke_color='#101028' , stroke_opacity=0 , fill_color='#101028' , fill_opacity=1 , width=config["frame_height"]*16.0/9.0 , height=config["frame_height"] , z_index=5)
        self.play(FadeIn(cover) , run_time=3 , func_rate=rush_into)
        self.clear()
        self.add(cover)
        self.play(Wait(run_time=1))
























class Part2a_2D(Scene):
    def construct(self):
        self.camera.background_color = '#101020'
        #create coordinate plane
        axes = ComplexPlane(
            x_range=(-3, 3, 1), y_range=(-3, 3, 1),
            x_length=6, y_length=6
        ).add_coordinates().set_opacity(0.25).scale(3).shift(0.3*DOWN+0.25*RIGHT)
        self.play(FadeIn(axes), run_time=2)
        self.wait(1)

        v1 = Arrow(axes.n2p(0+0j), axes.n2p(sqrt(3)/2 + 0.5j), color=RED, buff=0)
        v1_label = MathTex(r"\vec{v_1}", color=RED).next_to(v1, 0.5*UR)
        self.play(GrowArrow(v1), FadeIn(v1_label))
        self.wait(1)

        v2 = Arrow(axes.n2p(0+0j), axes.n2p(0.5 + 1j*sqrt(3)/2), color=BLUE, buff=0)
        v2_label = MathTex(r"\vec{v_2}", color=BLUE).next_to(v2, 0.5*UR)
        v3 = Arrow(axes.n2p(0+0j), axes.n2p(-sqrt(2)/2 - 1j*sqrt(2)/2), color=BLUE, buff=0)
        v3_label = MathTex(r"\vec{v_3}", color=BLUE).next_to(v3, 0.5*DL)
        self.play(GrowArrow(v2), GrowArrow(v3), FadeIn(v2_label), FadeIn(v3_label))
        self.wait(1)

        alpha_line = Angle(v1,v2, radius=0.5)
        beta_line = Angle(v1,v3, radius=0.3, other_angle=True)
        self.play(Create(alpha_line), Create(beta_line))

        alpha_label = MathTex(r"\alpha" , color=GREEN).next_to(alpha_line, 0.5*UR)
        beta_label = MathTex(r"\beta" , color=GREEN).next_to(beta_line, 0.5*DR)
        self.play(Write(alpha_label), Write(beta_label))
        self.wait(2)

        equation = MathTex(r"\alpha", r"<", r"\beta" , color=GREEN, font_size=144).move_to(axes.n2p(-1 + 0.9j))
        equation[1].set_fill(color=WHITE)
        self.play(Write(equation), run_time=2)

        conclusion = MathTex(r"\therefore D(", r"\vec{v_1}", r",", r"\vec{v_2}", r")", r"<", r"D(", r"\vec{v_1}", r",", r"\vec{v_3}", r")" , color=WHITE, font_size=48).move_to(axes.n2p(-1 + 0.4j))
        conclusion[1].set_fill(color=RED)
        conclusion[3].set_fill(color=BLUE)
        conclusion[7].set_fill(color=RED)
        conclusion[9].set_fill(color=BLUE)
        self.play(Write(conclusion), run_time=2)
        self.wait(5)

        self.play(Unwrite(conclusion))
        self.play(Unwrite(equation))
        self.wait(0.5)
        O = Dot(ORIGIN , radius=0)
        self.play(
            ReplacementTransform(v1,O),
            ReplacementTransform(v2,O),
            ReplacementTransform(v3,O),
            ReplacementTransform(v1_label,O),
            ReplacementTransform(v2_label,O),
            ReplacementTransform(v3_label,O),
            ReplacementTransform(alpha_label,O),
            ReplacementTransform(beta_label,O),
            ReplacementTransform(alpha_line,O),
            ReplacementTransform(beta_line,O),
            ReplacementTransform(equation,O),
            ReplacementTransform(conclusion,O),
        )
        self.wait(1)
        self.play(FadeOut(axes), run_time=1.5)
        self.wait(2.5)

















class Part1a_Projection(ThreeDScene):
    def construct(self):
        self.camera.background_color = '#101020'
        self.camera.light_source.move_to(OUT*20)
        self.set_camera_orientation(phi=2*PI/5, theta=PI/4)
        self.set_camera_orientation(zoom=0.35)
        axes = ThreeDAxes(
            x_range=(-12, 12, 2), y_range=(-12, 12, 2), z_range=(-8, 8, 2),
            x_length=24, y_length=24, z_length=16
        ).set_opacity(0.25)
        x_axis_label_rotation_angle , x_axis_label_rotation_axis = angle_axis_from_quaternion(quaternion_mult(quaternion_from_angle_axis(angle=PI, axis=OUT, axis_normalized=True) , quaternion_from_angle_axis(angle=PI/2, axis=RIGHT, axis_normalized=True)))
        y_axis_label_rotation_angle , y_axis_label_rotation_axis = angle_axis_from_quaternion(quaternion_mult(quaternion_from_angle_axis(angle=PI/2, axis=OUT, axis_normalized=True) , quaternion_from_angle_axis(angle=PI/2, axis=RIGHT, axis_normalized=True)))
        x_axis_label = axes.get_x_axis_label(MathTex("x" , color = RED).scale(3) , edge=RIGHT , buff=1).move_to([12.5,0,0]).rotate(x_axis_label_rotation_angle , axis=x_axis_label_rotation_axis)
        y_axis_label = axes.get_y_axis_label(MathTex("y" , color = GREEN).scale(3) , edge=UP , buff=1 , rotation=0).move_to([0,12.5,0]).rotate(y_axis_label_rotation_angle , axis=y_axis_label_rotation_axis)
        z_axis_label = axes.get_z_axis_label(MathTex("z" , color = BLUE).scale(3) , edge=OUT , buff=1 , rotation=PI/2).move_to([0,0,8.5]).rotate(PI*3/4 , axis=OUT)
        axes_labels = Group(x_axis_label , y_axis_label , z_axis_label)

        # define L1 and L2
        P1 = np.array([4,5,1])
        V1 = np.array([5,5,-4])
        P2 = np.array([4,-6,7])
        V2 = np.array([1,8,-3])
        N = get_unit_normal(V1,V2)
        
        # create L1 & L2
        L1 = Line3D([-4.75,-3.75,8] , [11,12,-4.6], thickness=0.05, color=PURE_RED)
        L2 = Line3D([11/3,-26/3,8] , [6.25,12,0.25], thickness=0.05, color=PURE_BLUE)

         # create Plane1
        Plane1 = Surface(
            lambda u, v: [u, v, 158/35 -(17/35)*u - (11/35)*v],
            u_range=[-12, 12],
            v_range=[-12, 12],
            resolution=32,
            fill_opacity=0.1,
            checkerboard_colors=[GREEN , PURE_GREEN],
            stroke_width=0
        )

        # create Plane2
        Plane2 = Surface(
            lambda u, v: [u, v, 247/35 -(17/35)*u - (11/35)*v],
            u_range=[-12, 12],
            v_range=[-12, 12],
            resolution=32,
            fill_opacity=0.1,
            checkerboard_colors=[GREEN , PURE_GREEN],
            stroke_width=0
        )

        # labels
        key_label = Tex(r"KEY" , color=WHITE).to_corner(UP + RIGHT)
        L1_label = MathTex(r"L_1" , color=PURE_RED).to_corner(UP + RIGHT).shift(DOWN*0.5)
        L2_label = MathTex(r"L_2" , color=PURE_BLUE).to_corner(UP + RIGHT).shift(DOWN*1)
        self.add_fixed_in_frame_mobjects(key_label, L1_label, L2_label)

        # fade in scene
        self.play(FadeIn(axes) , FadeIn(axes_labels) , FadeIn(L1) , FadeIn(L2) , FadeIn(Plane1) , FadeIn(Plane2) , FadeIn(key_label) , FadeIn(L1_label) , FadeIn(L2_label) , run_time=2)
        self.add(key_label, L1_label, L2_label)
        self.wait(0.5)
        self.move_camera(phi=9*PI/20, theta=-1*PI/50, zoom=0.35, run_time=1.5)
        self.wait(1)

        #create n
        n = Arrow3D(P2 - N*89/sqrt(1635) , P2 - N*89/sqrt(1635) + N, thickness=0.04, color=PURE_GREEN)

        n_label = MathTex(r"\hat{n}" , color=PURE_GREEN).to_corner(UP + RIGHT).shift(DOWN*1.5)
        self.add_fixed_in_frame_mobjects(n_label)
        self.play(GrowFromPoint(n, n.start), Write(n_label, reverse=True), run_time=1.5)
        self.add(n_label)
        self.wait(2)

        # create u
        u = Arrow3D(P1 , P2, thickness=0.05, color='#A000A0', z_index=2)

        u_label = MathTex(r"\vec{u}" , color='#A000A0').to_corner(UP + RIGHT).shift(DOWN*2)
        u_formula = MathTex(r"\vec{u}=\vec{p_2}-\vec{p_1}" , color='#A000A0').move_to(UP*1).rotate(-PI/6+0.03)
        self.add_fixed_in_frame_mobjects(u_label, u_formula)
        self.play(GrowFromPoint(u,P1), Write(u_label, reverse=True), Write(u_formula, reverse=True), run_time=1.5)
        self.add(u_label, u_formula)
        self.wait(2)
        self.play(Unwrite(u_formula, reverse=True))
        self.wait(2)

        #create u'
        uprime = Arrow3D(n.start , P2, thickness=0.05, color='#666666', z_index=3)

        uprime_label = MathTex(r"\vec{u'}" , color='#666666').to_corner(UP + RIGHT).shift(DOWN*2.5)
        uprime_formula = MathTex(r"\vec{u'}=\text{Proj}_{\hat{n}}\vec{u}" , color='#666666').move_to([-3,1.75,0]).rotate(-PI/60)
        self.add_fixed_in_frame_mobjects(uprime_label, uprime_formula)
        self.play(ReplacementTransform(u.copy(),uprime), Write(uprime_label, reverse=True), Write(uprime_formula, reverse=False), run_time=1.5)
        self.add(uprime_label, uprime_formula)
        self.wait(2)
        self.play(Unwrite(uprime_formula, reverse=False))
        self.wait(1)

        #move camera
        self.move_camera(phi=11*PI/20, theta=-1*PI/50, zoom=0.35, run_time=1.5)
        self.wait(1)

        #show D
        D_formula = MathTex(r"D=\| \vec{u'} \| \begin{cases} \\ \end{cases}").move_to([-3.75,2.95,0]).rotate(-PI/15)
        self.add_fixed_in_frame_mobjects(D_formula)
        self.play(Write(D_formula, reverse=True), run_time=1.5)
        self.add(D_formula)
        self.wait(2)
        self.play(Unwrite(D_formula, reverse=True))
        self.wait(1)

        #ending
        O = Dot3D(ORIGIN , radius=0)
        self.play(
            ReplacementTransform(key_label, O),
            ReplacementTransform(L1_label, O),
            ReplacementTransform(L2_label, O),
            ReplacementTransform(n_label, O),
            ReplacementTransform(u_label, O),
            ReplacementTransform(uprime_label, O), 
            ReplacementTransform(Plane1, O),
            ReplacementTransform(Plane2, O),
            ReplacementTransform(L1, O),
            ReplacementTransform(L2, O),
            ReplacementTransform(u, O),
            ReplacementTransform(n, O),
            ReplacementTransform(uprime, O),
        )
        self.wait(0.5)
        self.move_camera(phi=2*PI/5, theta=PI/4, zoom=0.35, run_time=1.5)
        self.wait(1)
        self.play(FadeOut(axes) , FadeOut(axes_labels) , run_time=1.5)
        self.wait(2.5)







class Part1a_Parallel_Planes(ThreeDScene):
    def construct(self):
        self.camera.background_color = '#101020'
        #self.camera.light_source.move_to(OUT*20)
        self.set_camera_orientation(phi=0, theta=-PI/2)
        self.set_camera_orientation(zoom=0.25)

        axes = ThreeDAxes(
            x_range=(-12, 12, 2), y_range=(-12, 12, 2), z_range=(-8, 8, 2),
            x_length=24, y_length=24, z_length=16
        ).set_opacity(0.25)
        x_axis_label = axes.get_x_axis_label(MathTex("x" , color = RED).scale(3) , edge=RIGHT , buff=1).move_to([12.5,0,0])
        y_axis_label = axes.get_y_axis_label(MathTex("y" , color = GREEN).scale(3) , edge=UP , buff=1 , rotation=0).move_to([0,12.5,0])
        z_axis_label = axes.get_z_axis_label(MathTex("z" , color = BLUE).scale(3) , edge=OUT , buff=1 , rotation=PI/2).move_to([0,0,8.5])
        axes_labels = Group(x_axis_label , y_axis_label , z_axis_label)
        self.play(FadeIn(axes) , FadeIn(axes_labels) , run_time=1.5)
        self.wait(0.5)

        # animate the move of the camera to properly see the axes
        self.move_camera(phi=2*PI/5, theta=-2*PI/5, zoom=0.35, run_time=1.5)
        self.begin_ambient_camera_rotation(rate=0.15)
        x_axis_label_rotation_angle , x_axis_label_rotation_axis = angle_axis_from_quaternion(quaternion_mult(quaternion_from_angle_axis(angle=PI, axis=OUT, axis_normalized=True) , quaternion_from_angle_axis(angle=PI/2, axis=RIGHT, axis_normalized=True)))
        y_axis_label_rotation_angle , y_axis_label_rotation_axis = angle_axis_from_quaternion(quaternion_mult(quaternion_from_angle_axis(angle=PI/2, axis=OUT, axis_normalized=True) , quaternion_from_angle_axis(angle=PI/2, axis=RIGHT, axis_normalized=True)))
        self.play(Rotate(axes_labels[0] , x_axis_label_rotation_angle , axis=x_axis_label_rotation_axis),
                  Rotate(axes_labels[1] , y_axis_label_rotation_angle , axis=y_axis_label_rotation_axis),
                  Rotate(axes_labels[2] , PI*3/4 , axis=OUT) ,
                  run_time=1)
        self.wait(1.5)

        # define L1 and L2
        P1 = np.array([4,5,1])
        V1 = np.array([5,5,-4])
        P2 = np.array([4,-6,7])
        V2 = np.array([1,8,-3])

        # create L1
        L1 = Line3D([-4.75,-3.75,8] , [11,12,-4.6], thickness=0.02, color=PURE_RED)
        L1_label = MathTex(r"L_1:\begin{pmatrix}x\\y\\z\end{pmatrix}=\begin{pmatrix}4\\5\\1\end{pmatrix}+t\begin{pmatrix}5\\5\\-4\end{pmatrix}").rotate(angle=PI/2, axis=RIGHT).rotate(angle=PI/4, axis=OUT).move_to([5,5,5]).scale(2)
        self.play(Create(L1), FadeIn(L1_label))
        self.wait(3)
        self.play(FadeOut(L1_label))
        self.wait(1)

        # create L2
        L2 = Line3D([11/3,-26/3,8] , [6.25,12,0.25], thickness=0.02, color=PURE_BLUE)
        L2_label = MathTex(r"L_2:\begin{pmatrix}x\\y\\z\end{pmatrix}=\begin{pmatrix}4\\-6\\7\end{pmatrix}+s\begin{pmatrix}1\\8\\-3\end{pmatrix}").rotate(angle=PI/2, axis=RIGHT).rotate(angle=PI/2, axis=OUT).move_to([3,6,5]).scale(2)
        self.play(Create(L2), FadeIn(L2_label))
        self.wait(3)
        self.play(FadeOut(L2_label))
        self.wait(1)

        # get V1 and V2
        V1_vector = Arrow3D(P1, P1+V1 , color=PURE_RED)
        V2_vector = Arrow3D(P2, P2+V2 , color=PURE_BLUE)
        V1_label = MathTex(r"\vec{v_1}" , color=PURE_RED).rotate(angle=PI/2, axis=RIGHT).rotate(angle=PI, axis=OUT).move_to(V1+[1,1,0]).scale(3)
        V2_label = MathTex(r"\vec{v_2}" , color=PURE_BLUE).rotate(angle=PI/2, axis=RIGHT).rotate(angle=PI, axis=OUT).move_to(V2+[-1,1,0]).scale(3)
        self.play(FadeIn(V1_vector), FadeIn(V2_vector))
        self.play(ApplyMethod(V1_vector.shift, -P1), FadeIn(V1_label))
        self.wait(0.5)
        self.play(ApplyMethod(V2_vector.shift, -P2), FadeIn(V2_label))
        self.wait(1.5)
        self.play(FadeOut(V1_label), FadeOut(V2_label))
        self.wait(1)

        # create Plane0
        Plane0 = Surface(
            lambda u, v: [u, v, -(17/35)*u - (11/35)*v],
            u_range=[-12, 12],
            v_range=[-12, 12],
            resolution=64,
            fill_opacity=0.5,
            checkerboard_colors=[GREEN , PURE_GREEN],
        )
        self.play(GrowFromPoint(Plane0, ORIGIN))
        self.wait(3)

        # create Plane1
        Plane1 = Surface(
            lambda u, v: [u, v, 158/35 -(17/35)*u - (11/35)*v],
            u_range=[-12, 12],
            v_range=[-12, 12],
            resolution=64,
            fill_opacity=0.5,
            checkerboard_colors=[GREEN , PURE_GREEN],
        )
        self.play(ReplacementTransform(Plane0,Plane1))
        self.wait(2)

        # create Plane2
        Plane2 = Surface(
            lambda u, v: [u, v, 247/35 -(17/35)*u - (11/35)*v],
            u_range=[-12, 12],
            v_range=[-12, 12],
            resolution=64,
            fill_opacity=0.5,
            checkerboard_colors=[GREEN , PURE_GREEN],
        )
        self.play(ReplacementTransform(Plane1,Plane2))
        self.wait(2)

        #closing
        Plane0 = Surface(
            lambda u, v: [u, v, -(17/35)*u - (11/35)*v],
            u_range=[-12, 12],
            v_range=[-12, 12],
            resolution=64,
            fill_opacity=0.5,
            checkerboard_colors=[GREEN , PURE_GREEN],
        )
        self.play(ReplacementTransform(Plane2,Plane0))
        self.wait(5)

        O = Dot3D(ORIGIN , radius=0)
        self.play(
            ReplacementTransform(Plane0,O),
            ReplacementTransform(L1,O),
            ReplacementTransform(L2,O),
            ReplacementTransform(V1_vector,O),
            ReplacementTransform(V2_vector,O)
        )
        self.wait(0.5)
        self.play(Rotate(axes_labels[0] , -x_axis_label_rotation_angle , axis=x_axis_label_rotation_axis),
                  Rotate(axes_labels[1] , -y_axis_label_rotation_angle , axis=y_axis_label_rotation_axis),
                  Rotate(axes_labels[2] , -PI*3/4 , axis=OUT) ,
                  run_time=1)
        self.stop_ambient_camera_rotation()
        self.move_camera(phi=0, theta=-PI/2, zoom=0.25, run_time=1.5)
        self.wait(1)
        self.play(FadeOut(axes) , FadeOut(axes_labels) , run_time=1.5)
        self.wait(2.5)