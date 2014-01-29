
//        double H = 0;
//        double G = 0;
//        Node integ = integral(p, a, b);
//        G = integ.y; // Analytic solution
//        H = integ.x; // Analytic solution


Node Bem::integral(Node v, Node pa, Node pb)
{
    // center
    Node center = (pa + pb) / 2;

    // shift of center of the edge on the position [0, 0]
    Node at = pa - center;
    Node bt = pb - center;

    // shift of reference point
    Node vt = v - center;

    // rotation
    double phi = - atan2(bt.y, bt.x);
    vt.rotate(phi);
    at.rotate(phi);
    bt.rotate(phi);

    double H;
    double G;

    if(v.distanceOf(center) < EPS_ZERO)
    {
        // singular point
        H = 0;
        double m = abs(bt.x);
        G = -1 / ( M_PI)  * (m * log(m) - m);
        return Node(H, G);
    } else
    {

        double x = 0;
        double y = 0;
        double a = (vt.x - at.x);
        double b = (vt.x - bt.x);

        if(abs(vt.y) < EPS_ZERO)
        {

            if((abs(vt.x - at.x) < EPS_ZERO) && (abs(vt.x - bt.x) < EPS_ZERO))
            {
                assert(0);
            }

            if(abs(b) < EPS_ZERO)
            {
                x = (a > 0) ? M_PI_2 : - M_PI_2;
                y = M_PI_4;
                G =  - ( a * (-2 + log(a * a)))/(4 * M_PI);
            } else
                if(abs(a) < EPS_ZERO)
                {
                    y = (b > 0) ? M_PI_2 : - M_PI_2;
                    x = M_PI_4;
                    G = (b * (-2 + log(b * b)))/(4 * M_PI);
                }
                else
                {
                    x = (a > 0) ? M_PI_2 : - M_PI_2;
                    y = (b > 0) ? M_PI_2 : - M_PI_2;
                    G = (b * (-2 + log(b * b)))/(4 * M_PI) - (a * (-2 + log(a * a)))/(4 * M_PI);
                }
        }
        else
        {
            x = atan( a / vt.y);
            y = atan( b / vt.y);
            G = (2 * vt.y * y + b * (-2 + log(vt.y * vt.y + b * b)))/(4 * M_PI) -
                    (2 * vt.y * x + a * (-2 + log(vt.y * vt.y + a * a)))/(4 * M_PI);
        }
        H = (y - x)  / (2 * M_PI);
    }

    return Node(H, G);
}
