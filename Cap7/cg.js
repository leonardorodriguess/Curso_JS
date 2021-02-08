class Point2D {
    constructor(x, y){
        this.x = x;
        this.y = y;
    }

    Add(vector){
        return new Point2D(vector.x + this.x, vector.y + this.y);
    }

    Sub(point){
        return new Vector2D(this.x - point.x, this.y - point.y);
    }
}

class Point3D {
    constructor(x, y, z, w = 1.0){
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
    }

    Add(vector){
        return new Point3D(vector.x + this.x, vector.y + this.y, vector.z + this.z);
    }

    Sub(point){
        return new Vector3D(this.x - point.x, this.y - point.y, this.z - point.z);
    }

    ToCart() {
        return new Point3D(this.x/this.w, this.y/this.w, this.z/this.w, 1.0);
    }
}

class Vector2D {
    constructor(dx, dy) {
        this.dx = dx;
        this.dy = dy;
    }

    normalize(){
        let n = Math.sqrt(this.dx * this.dx + this.dy * this.dy);
        this.dx /= n;
        this.dy /= n;
        return this;
    }
}

class Vector3D {
    constructor(dx, dy, dz) {
        this.x = dx;
        this.y = dy;
        this.z = dz;
        this.w = 0;
    }

    normalize(){
        let n = Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
        this.x /= n;
        this.y /= n;
        this.z /= n;
        return this;
    }

    scale(k){
        return new Vector3D(this.x * k, this.y * k, this.z * k);
    }

    dot(v){
        return this.x * v.x + this.y * v.y + this.z * v.z; 
    }

    cross(v){
        return new Vector3D(this.y * v.z - this.z * v.y, this.z * v.x - this.x * v.z, this.x * v.y - this.y * v.x);
    }
}

class Camera {
    constructor(eye, center, up) {
        this.eye = eye;
        this.w = center.Sub(eye).normalize();
        this.u = up.cross(this.w).normalize();
        this.v = this.w.cross(this.u).normalize();
    }
}

class Matrix33 {
    constructor(a00=1, a01=0, a02=0, a10=0, a11=1, a12=0, a20=0, a21=0, a22=1) {
        this.a00 = a00;
        this.a01 = a01;
        this.a02 = a02;
        this.a10 = a10;
        this.a11 = a11;
        this.a12 = a12;
        this.a20 = a20;
        this.a21 = a21;
        this.a22 = a22;
    }

    mult(m){
        let r = new Matrix33();
        r.a00 = this.a00 * m.a00 + this.a01 * m.a10 + this.a02 * m.a20;
        r.a01 = this.a00 * m.a01 + this.a01 * m.a11 + this.a02 * m.a21;
        r.a02 = this.a00 * m.a02 + this.a01 * m.a12 + this.a02 * m.a22;
    
        r.a10 = this.a10 * m.a10 + this.a11 * m.a10 + this.a12 * m.a20;
        r.a11 = this.a10 * m.a01 + this.a11 * m.a11 + this.a12 * m.a21;
        r.a12 = this.a10 * m.a02 + this.a11 * m.a12 + this.a12 * m.a22;

        r.a20 = this.a20 * m.a00 + this.a21 * m.a10 + this.a22 * m.a20;
        r.a21 = this.a20 * m.a01 + this.a21 * m.a11 + this.a22 * m.a21;
        r.a22 = this.a20 * m.a02 + this.a21 * m.a12 + this.a22 * m.a22;
        return r;
    }

    det() {
        let a = this.a20 * this.a01 * this.a12 + this.a00 * this.a11 * this.a22 + this.a10 * this.a21 * this.a02; 
        let b = this.a00 * this.a21 * this.a12 + this.a20 * this.a11 * this.a02 + this.a10 * this.a01 * this.a22;
        return a - b;
    }
}


class Matrix44 {
    constructor(a00=1, a01=0, a02=0, a03=0, a10=0, a11=1, a12=0, a13=0, a20=0, a21=0, a22=1, a23=0, a30=0, a31=0, a32=0, a33=0) {
        this.a00 = a00;
        this.a01 = a01;
        this.a02 = a02;
        this.a03 = a03;
        this.a10 = a10;
        this.a11 = a11;
        this.a12 = a12;
        this.a13 = a13;
        this.a20 = a20;
        this.a21 = a21;
        this.a22 = a22;
        this.a23 = a23;
        this.a30 = a30;
        this.a31 = a31;
        this.a32 = a32;
        this.a33 = a33;
    }

    mult(m){
        let r = new Matrix44();
        r.a00 = this.a00 * m.a00 + this.a01 * m.a10 + this.a02 * m.a20 + this.a03 * m.a30;
        r.a01 = this.a00 * m.a01 + this.a01 * m.a11 + this.a02 * m.a21 + this.a03 * m.a31;
        r.a02 = this.a00 * m.a02 + this.a01 * m.a12 + this.a02 * m.a22 + this.a03 * m.a32;
        r.a03 = this.a00 * m.a03 + this.a01 * m.a13 + this.a02 * m.a23 + this.a03 * m.a33;
    

        r.a10 = this.a10 * m.a00 + this.a11 * m.a10 + this.a12 * m.a20 + this.a13 * m.a30;
        r.a11 = this.a10 * m.a01 + this.a11 * m.a11 + this.a12 * m.a21 + this.a13 * m.a31;
        r.a12 = this.a10 * m.a02 + this.a11 * m.a12 + this.a12 * m.a22 + this.a13 * m.a32;
        r.a13 = this.a10 * m.a03 + this.a11 * m.a13 + this.a12 * m.a23 + this.a13 * m.a33;

        r.a20 = this.a20 * m.a00 + this.a21 * m.a10 + this.a22 * m.a20 + this.a23 * m.a30;
        r.a21 = this.a20 * m.a01 + this.a21 * m.a11 + this.a22 * m.a21 + this.a23 * m.a31;
        r.a22 = this.a20 * m.a02 + this.a21 * m.a12 + this.a22 * m.a22 + this.a23 * m.a32;
        r.a23 = this.a20 * m.a03 + this.a21 * m.a13 + this.a22 * m.a23 + this.a23 * m.a33;

        r.a30 = this.a30 * m.a00 + this.a31 * m.a10 + this.a32 * m.a20 + this.a33 * m.a30;
        r.a31 = this.a30 * m.a01 + this.a31 * m.a11 + this.a32 * m.a21 + this.a33 * m.a31;
        r.a32 = this.a30 * m.a02 + this.a31 * m.a12 + this.a32 * m.a22 + this.a33 * m.a32;
        r.a33 = this.a30 * m.a03 + this.a31 * m.a13 + this.a32 * m.a23 + this.a33 * m.a33;

        return r;
    }

    transform(p, isPoint=true) {
        let xr = this.a00 * p.x + this.a01 * p.y + this.a02 * p.z + this.a03 * p.w;
        let yr = this.a10 * p.x + this.a11 * p.y + this.a12 * p.z + this.a13 * p.w;
        let zr = this.a20 * p.x + this.a21 * p.y + this.a22 * p.z + this.a23 * p.w;
        let wr = this.a30 * p.x + this.a31 * p.y + this.a32 * p.z + this.a33 * p.w;
        if (isPoint){
            return new Point3D(xr, yr, zr, wr);
        } else {
            return new Vector3D(xr, yr, zr, wr);
        }
    }

    toString(d=2) {
        return this.a00.toFixed(d) + ", " + this.a01.toFixed(d) + ", " + this.a02.toFixed(d) + ", " + this.a03.toFixed(d) + "\n" +
               this.a10.toFixed(d) + ", " + this.a11.toFixed(d) + ", " + this.a12.toFixed(d) + ", " + this.a13.toFixed(d) + "\n" +
               this.a20.toFixed(d) + ", " + this.a21.toFixed(d) + ", " + this.a22.toFixed(d) + ", " + this.a23.toFixed(d) + "\n" +
               this.a30.toFixed(d) + ", " + this.a31.toFixed(d) + ", " + this.a32.toFixed(d) + ", " + this.a33.toFixed(d) + "\n";
    }
}


class Line3D {
    constructor(begin, end, r=1, g=0, b=0){
        this.begin = begin;
        this.end = end;
        this.r = r;
        this.g = g;
        this.b = b;
    }
}

class WiredQuad {
    constructor(a, b, c, d, r=1, g=0, cb=0){
        this.lines = [new Line3D(a, b, r, g, cb), 
                      new Line3D(b, c, r, g, cb), 
                      new Line3D(c, d, r, g, cb),
                      new Line3D(d, a, r, g, cb)]
        this.r = r;
        this.b = b;
        this.g = g;
    }
}

class WiredTriangle {
    constructor(a, b, c, r=0, g=0, cb=0) {
        this.lines = [new Line3D(a, b, r, g, cb), new Line3D(b, c, r, g, cb), new Line3D(c, a, r, g, cb)];
    }
}


class WiredSphere {
    constructor(pos, radio, segments, rings){
        this.radio = radio;
        this.position = pos;
        this.faces = [];
        let theta = 0;
        let deltaTheta = Math.PI/rings;
        let phi = -Math.PI;
        let deltaPhi = (2 * Math.PI)/segments;
        
        let ring0 = [new Point3D(pos.x, pos.y, radio+pos.z)];
        for (let t = theta+deltaTheta; t < Math.PI; t += deltaTheta) {
            let ringN = []
            let z = radio * Math.cos(t) + pos.z;
            for (let p = phi; p < Math.PI; p += deltaPhi) {
                let x = radio * Math.cos(p) * Math.sin(t) + pos.x;
                let y = radio * Math.sin(p) * Math.sin(t) + pos.y;
                ringN.push(new Point3D(x, y, z));
            }
            let x = radio * Math.cos(Math.PI) * Math.sin(t) + pos.x;
            let y = radio * Math.sin(Math.PI) * Math.sin(t) + pos.y;
            ringN.push(new Point3D(x, y, z));

            if (ring0.length == 1) {
                for (let j = 0; j < ringN.length; j++) {
                    let a = j;
                    let b = j + 1;
                    if (b >= ringN.length) {
                        b = 0;
                    }
                    this.faces.push(new WiredTriangle(ring0[0], ringN[a], ringN[b]));
                }
            } else if (ring0.length > 1) {
                for (let j = 0; j < ringN.length; j++) {
                    let a = j;
                    let b = j + 1;
                    if (b >= ringN.length) {
                        b = 0;
                    }
                    this.faces.push(new WiredQuad(ring0[a], ringN[a], ringN[b], ring0[b]));
                }
            }
            ring0 = ringN;
        }
        let ringN = ring0
        ring0 = [new Point3D(pos.x, pos.y, -radio+pos.z)];
        for (let j = 0; j < ringN.length; j++) {
            let a = j;
            let b = j + 1;
            if (b >= ringN.length) {
                b = 0;
            }
            this.faces.push(new WiredTriangle(ring0[0], ringN[a], ringN[b]));
        }

        this.lines = []
        for (let i = 0; i < this.faces.length; i++) {
            for (let j = 0; j < this.faces[i].lines.length; j++) {
                this.lines.push(this.faces[i].lines[j]);
            }
        }
    }
}

class WiredCube {
    constructor(origin, size){
        this.origin = origin;
        this.size = size;

        let v1 = new Point3D(origin.x - size * 0.5, origin.y - size * 0.5, origin.z + size * 0.5);
        let v2 = new Point3D(origin.x + size * 0.5, origin.y - size * 0.5, origin.z + size * 0.5);
        let v3 = new Point3D(origin.x + size * 0.5, origin.y + size * 0.5, origin.z + size * 0.5);
        let v4 = new Point3D(origin.x - size * 0.5, origin.y + size * 0.5, origin.z + size * 0.5);

        let v5 = new Point3D(origin.x - size * 0.5, origin.y - size * 0.5, origin.z - size * 0.5);
        let v6 = new Point3D(origin.x + size * 0.5, origin.y - size * 0.5, origin.z - size * 0.5);
        let v7 = new Point3D(origin.x + size * 0.5, origin.y + size * 0.5, origin.z - size * 0.5);
        let v8 = new Point3D(origin.x - size * 0.5, origin.y + size * 0.5, origin.z - size * 0.5);
    
        this.faces = [new WiredQuad(v1, v2, v3, v4, 0, 255, 127), //front 
                      new WiredQuad(v5, v8, v7, v6, 255, 255, 255), //back
                      new WiredQuad(v1, v5, v6, v2, 255, 100, 0), //bottom
                      new WiredQuad(v4, v3, v7, v8, 190, 190, 190), //top
                      new WiredQuad(v1, v4, v8, v5, 0, 100, 255), //left
                      new WiredQuad(v2, v6, v7, v3, 255, 0, 0), //right
                    ]
        this.lines = []
        for (let i = 0; i < this.faces.length; i++) {
            this.lines.push(this.faces[i].lines[0]);
            this.lines.push(this.faces[i].lines[1]);
            this.lines.push(this.faces[i].lines[2]);
            this.lines.push(this.faces[i].lines[3]);
        }
    }
}

//implementa a visualização em wireframe
function render(context, seg) {
    let nx = 600; //eh a quantidade de pixels na horizontal.
    let ny = 600; //eh a quantidade de pixels na vertical.
    let n = -100; //a posicao do plano mais proximo da camera.
    let f = -200; //a posicao do plano mais distante da camera.
    let l = -100; //valor minimo de x
    let r = 100; //valor maximo de x
    let b = -100; //valor minimo de y
    let t = 100; //valor maximo de y
    let camx = 10, camy = 10, camz = 10; //posicao da camera em coordenadas do mundo

    context.save()
    context.clearRect(0,0, 600,600)

    let camera = new Camera(new Point3D(seg*Math.sin(seg*Math.PI/180),seg*Math.cos(seg*Math.PI/180), 10), new Point3D(0, 0, 1), new Vector3D(1, 0, 1))
    let u = camera.u
    let v = camera.v
    let w = camera.w
    let eye = camera.eye
    
    let mvp =  new Matrix44(nx/2.0, 0, 0, (nx-1)/2.0,       0, ny/2.0, 0, (ny-1)/2.0,    0, 0, 1, 0,    0, 0, 0, 1);
    let mortho = new Matrix44(2/(r-l), 0, 0, -(r+l)/(r-l),  0, 2/(t-b), 0, -(t+b)/(t-b),  0, 0, 2/(n-f), -(n+f)/(n-f),   0, 0, 0, 1.0);
    let camrot = new Matrix44(u.x, u.y, u.z, 0,   v.x, v.y, v.z, 0,   w.x, w.y, w.z, 0,   0, 0, 0, 1)
    let campos = new Matrix44(1, 0, 0, -eye.x,   0, 1, 0, -eye.y,   0, 0,  1, -eye.z,   0, 0, 0, 1)
    let mcam = camrot.mult(campos)
    
    let mview = mvp.mult(mortho).mult(mcam)
    
    console.log(mview.toString());
    

    let objects = [new WiredCube(new Point3D(0, 0, 10), 50)];

    let i = 0
    while (i < objects.length){
        let obj = objects[i]
        for (j = 0; j < obj.lines.length; j++) {
            line = obj.lines[j]
            let a = mview.transform(line.begin);
            //a = a.ToCart();
            let b = mview.transform(line.end);
            //b = b.ToCart();
            //console.log("(" + a.x + ", " + a.y + ", " + a.z + ") --- (" + b.x + ", " + b.y + ", " + b.z + ")");
            DrawLine(context, a.x, a.y, b.x, b.y, line.r, line.g, line.b, 255, 0,
                    nx-1, 0, ny-1, nx, ny)
            /*let bx = b.x 
            let ay = a.y
            let by = b.y
            for(let ax = a.x; ax < a.y; ax += 0.5){
                bx += 0.5
                ay += 0.5
                by += 0.5
                DrawLine(context, bx, a.y, bx, b.y, line.x, line.g, line.b, 255, 0,
                    nx-1, 0, ny-1, nx, ny, a.z, b.z)

            }*/

                
                
            

            
        }
        i++;
    }
}


function XToScr(X, MINXW=-100, MAXXW=100, MINXS=0, MAXXS=300)
{
	let XN = (X-MINXW)/(MAXXW - MINXW);
    return Math.trunc(MINXS + XN * (MAXXS - MINXS) + 0.5);
}

function YToScr(Y, MINYW=-100, MAXYW=100, MINYS=0, MAXYS=150)
{
	let YN = (Y-MINYW)/(MAXYW - MINYW);
    return Math.trunc(MAXYS - YN * (MAXYS - MINYS) + 0.5);
}

function XToWrld(X, MINXW=-100, MAXXW=100, MINXS=0, MAXXS=300) {
	let XN = (X-MINXS)/(MAXXS - MINXS);
    return MINXW + XN * (MAXXW - MINXW);
}

function YToWrld(Y, MINYW=-100, MAXYW=100, MINYS=0, MAXYS=150)
{
	let YN = (Y-MINYS)/(MAXYS - MINYS);
    return Math.trunc(MAXYW - YN * (MAXYW - MINYW));
}

function DrawLine(context, x1w, y1w,  x2w, y2w, r=0, g=0, b=0, a=255, xmin=-1, xmax=1, ymin=-1, ymax=1, nx = 600, ny = 600, z1w = 0, z2w = 0) {
    let bkp = context.strokeStyle;
    context.strokeStyle = "rgba(" + r + "," + g + "," + b + "," + (a/255.0) + ")";
    context.beginPath(); 
    // Staring point (10,45)
     context.moveTo(XToScr(x1w, xmin, xmax, 0, nx, 0, ny), YToScr(y1w, ymin, ymax, 0, nx, 0, ny));
    // End point (180,47)
    context.lineTo(XToScr(x2w, xmin, xmax, 0, nx, 0, ny), YToScr(y2w, ymin, ymax, 0, nx, 0, ny));
    // Make the line visible
    context.stroke();
    context.strokeStyle = bkp;
}

function DrawMarker(targetContext, xw, yw, r=0, g=0, b=0, a=255, size = 10)
{
    let bkp = targetContext.fillStyle;
    targetContext.fillStyle = "rgba(" + r + "," + g + "," + b + "," + (a/255.0) + ")";
    targetContext.fillRect(XToScr(xw-size/2.0), YToScr(yw+size/2.0), size, size);
    targetContext.fillStyle = bkp;
}

function SetPixelAtWorld(targetContext, xw, yw, r=0, g=0, b=0, a=255)
{
    let bkp = targetContext.fillStyle;
    targetContext.fillStyle = "rgba(" + r + "," + g + "," + b + "," + (a/255.0) + ")";
    targetContext.fillRect(XToScr(xw), YToScr(yw), 1, 1);
    targetContext.fillStyle = bkp;
}
