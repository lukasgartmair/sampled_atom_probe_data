import bpy
import mathutils

def createCurveObject():
    # Create curve and object
    cu = bpy.data.curves.new('MyCurve', 'CURVE')
    ob = bpy.data.objects.new('MyCurveObject', cu)
    bpy.context.scene.objects.link(ob)

#co <Vector (-0.3138, -0.0000, 0.0024)>
#handle_left <Vector (-0.3566, -0.2274, -0.6022)>
#handle_right <Vector (-0.2109, 0.5462, 1.4543)>
#co <Vector (5.8808, -0.0000, 1.7996)>
#handle_left <Vector (4.5356, 0.0000, 1.3798)>
#handle_right <Vector (6.2238, -0.0000, 1.9067)>

    P0 = (-0.3138, -0.0000, 0.0024)
    P1 = (-0.3566, -0.2274, -0.6022)
    P2 = (-0.2109, 0.5462, 1.4543)

    P3 = (5.8808, -0.0000, 1.7996)
    P4 = (4.5356, 0.0000, 1.3798)
    P5 = (6.2238, -0.0000, 1.9067)

    # Bezier coordinates
    beziers = [
                ((-0.3138, -0.0000, 0.0024), (-0.3566, -0.2274, -0.6022), (-0.2109, 0.5462, 1.4543))
                , ((5.8808, -0.0000, 1.7996), (4.5356, 0.0000, 1.3798), (6.2238, -0.0000, 1.9067))
              ]
 
    # Create spline and set Bezier control points
    spline = cu.splines.new('BEZIER')
    nPointsU = len(beziers)
    print(nPointsU)
    spline.bezier_points.add(nPointsU-1)
    for n in range(nPointsU):
        
        print(str(n) + 'n')
        bpt = spline.bezier_points[n]

        bpt.co = mathutils.Vector(beziers[n][0])
        bpt.handle_left = mathutils.Vector(beziers[n][1])
        bpt.handle_right = mathutils.Vector(beziers[n][2])
            
    return ob

curveob = createCurveObject()
curveob.select = True

bpy.context.scene.objects.active = curveob
obj = bpy.context.active_object

if obj.type == 'CURVE':
    for subcurve in obj.data.splines:
        curvetype = subcurve.type
        print('curve type:', curvetype)

        if curvetype == 'BEZIER':
            print("curve is closed:", subcurve.use_cyclic_u)

            # print(dir(subcurve))
            for bezpoint in subcurve.bezier_points:
                """
                'handle_left_type',      # kind of handles 
                'handle_right_type',     # 
                'hide',                  # is it hidden?
                'radius',                # what's the radius
                'select_control_point',  # is it selected?
                'select_left_handle',    #
                'select_right_handle',   #
                'tilt'                   # investigate :)
                # use     print(dir(bezpoint))  to see all
                """
                print('co', bezpoint.co)
                print('handle_left', bezpoint.handle_left)
                print('handle_right', bezpoint.handle_right)
                
bpy.context.scene.objects.active = curveob
act_obj = bpy.context.active_object
screw = act_obj.modifiers.new('screw', 'SCREW')
screw.axis = 'X'

bpy.ops.object.convert(target='MESH', keep_original=False)
