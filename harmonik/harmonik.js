function submit() {
 let    xmin    = parseFloat(document.getElementById('xmin').value),
        xmax    = parseFloat(document.getElementById('xmax').value),
        dx      = parseFloat(document.getElementById('dx').value),
        initE   = parseFloat(document.getElementById('initE').value),

        v0      = parseFloat(document.getElementById('v0').value),
        om      = parseFloat(document.getElementById('omega').value);
        ef      = document.getElementById('ef');

let x       = [],
    vpot    = [],
    idx     = [],
    efunction = [],
    vtemp   = [],
    iter    = 0,
    nloop   = 10,

    leb     = xmax - xmin,
    ngrid   = leb/dx;

    
    //console.log(initE + 1e-16)
    singleEigen(xmin,dx,ngrid,v0,vpot,x,idx,efunction,iter,nloop,initE,leb,om, vtemp,potensial_function = harmonicPotential)
    //console.log(new Float32Array(initE))

    

    // efunction.forEach(function(item) {
    //   var listItem = document.createElement('li');
    //   listItem.className = "nostyle";
    //   listItem.innerText = item;
    //   document.getElementById("ef").appendChild(listItem);
    //  });

     
}

// xmin = xmin.value;
// console.log(xmin,xmax,dx,initE,v0,omega);


    
function harmonicPotential(xmin,leb,om,x,v0 = 0)
{
    //$xp = [];
    xp = x - xmin - leb / 2;
    console.log(xp, x, xmin, leb)
    pot = (0.5 * om * xp * xp) + v0;
    return pot;
}

function singleEigen(xmin,dx,ngrid,v0,vpot,x,idx,efunction,iter,nloop,initE,leb, om, vtemp,potensial_function = harmonicPotential)
{
    for (let i = 0; i < ngrid; i++) {
        x[i] = xmin + ((i)  * dx);
        //console.log(x[i]);

        vpot[i] = harmonicPotential(xmin,leb,om,x[i], v0);
        idx[i] = i;
        //console.log(vpot[i])
    }
    // console.log(xmin,xmax,ngrid)
    //console.log(x);
    eigenState(x,ngrid,efunction,vpot,initE,iter,nloop,dx,vtemp);
}

function spaceship(val1, val2) {
    if ((val1 === null || val2 === null) || (typeof val1 != typeof val2)) {
        return null;
    }
    if (typeof val1 === 'string') {
        return (val1).localeCompare(val2);
    }
    else {
        if (val1 > val2) { return 1 }
        else if (val1 < val2) { return -1 }
        return 0;
    }
}

function eigenState(x,ngrid,efunction,vpot,initE,iter,nloop,dx,vtemp)
{
    //initE = initE
    //console.log(initE)
    initE = initE + 1e-16;
    
    //console.log(initE)
    x_tart = x[0];
    x_end  = x[ngrid - 1];
    if (x_tart < 0) {
        for (let i = 0; i < ngrid; i++) {
            efunction[i] = (Math.sin(x[i]) + Math.cos(x[i]));
        }

    } else {
        for (let i = 0; i < ngrid;i++) {
            efunction[i] = 1 + x[i] / x_end;
        }
    }

    vharm = vpot[ngrid - 1];
    estart = initE;
    eps = 1e-5;
    if (estart > 0) eps = 1e-5;
    if (vharm > 1) eps = 1e-8;
    energy = initE;
    isig = 1;
    if (initE < 0) isig = -1;
    do {
        //console.log(initE);
        iter += iter;
        invers(dx,vtemp,vpot,efunction,ngrid,nloop);
        sum = 0;
        for (let i = 0; i < ngrid; i++) {
            sum = sum + efunction[i] * efunction[i];
        }
        sum = Math.sqrt(sum * dx);
        for (i = 1; i < ngrid; i++) {
            efunction[i] = efunction[i] / sum;
        }
        if (x_tart > 0) hamilton(dx,vtemp,vpot,efunction,ngrid);
        else hamilton5p(dx,vtemp,vpot,efunction,ngrid);
        //var_dump(global vtemp);
        sum = 0;
        //console.log(vtemp);
        for (let i = 1; i < ngrid - 1; i++) {
            sum = sum + vtemp[i] * efunction[i];
        }
        energy_new = sum * dx;
         
        delta = Math.abs((energy_new - energy) / energy);

        energyN = energy_new;
        if (iter < 1) energyN = initE;
        energy = energyN;
    } while (delta <= eps);

    initE = energy_new;
   
    // vtemp = [];
    let rsign = 1;
    if (x < 1e-8) rsign = (spaceship(efunction[1],0));
    ampl_max = 1;
    if (vharm <= 1) {
        if (estart > 0 && initE > 0) {
            ampl_max = 0;
            inode = 0;
            wb = efunction[ngrid - 1];
            i = ngrid;
            while (inode < 10) {
                i = i - 1;
                wf = efunction[i];
                if (i < 10) break;
                if (wf * wb < 0) inode = inode + 1;
                abswf = Math.abs(wf);
                if (abswf > ampl_max) ampl_max = abswf;
                wb = wf;
            }
        }
    }
    ampl_max = ampl_max * rsign;
    for (let i = 0; i < ngrid; i++) {
        efunction[i] = efunction[i] / ampl_max;
    }
    //var_dump(efunction);
    console.log(initE);
    // for (let i = 0; i < efunction.length; i++) {
    // return (
    //     cob.innerHTML = "hihia"
    // )
    

    // //     console.log(efunction[i]);
    // }
      //ef.innerHTML = efunction;

    console.log(efunction)

    Enew.innerHTML = initE;

    efunction.forEach(function(item) {
        var listItem = document.createElement('li');
        listItem.className = "nostyle";
        listItem.innerText = item;
        document.getElementById("ef").appendChild(listItem);
     });
}

function hamilton(dx,vtemp,vpot,efunction,ngrid)
{
    let dtr = 1 / (dx * dx) / 2;
    let a = -dtr;
    let b = 2 * dtr;
    let c = -dtr;
    vtemp[0] = (b + vpot[0]) * efunction[0] + c * efunction[1];
    for (let i = 1; i < ngrid; i++) {
        vtemp[i] = (b + vpot[i]) * efunction[i] + c * (efunction[i - 1] + efunction[i + 1]);
    }
    vtemp[ngrid] = (b + vpot[ngrid]) * efunction[ngrid] + a * efunction[ngrid - 1];
    console.log(vtemp);
}
function hamilton5p(dx,vtemp,vpot,efunction,ngrid)
{
    let dtr = 1 / (dx * dx) / 24;
    let a = dtr;
    let b = 30 * dtr;
    let c = -16 * dtr;
    let i = 0;
    vtemp[i] = (b + vpot[i]) * efunction[i]
        + c * efunction[i + 1]
        + a * efunction[i + 2];
    i = 1;
    vtemp[i] = (b + vpot[i]) * efunction[i]
        + c * (efunction[i - 1] + efunction[i + 1])
        + a * efunction[i + 2];
    for (let i = 2; i < ngrid - 3; i++) {
        vtemp[i] = (b + vpot[i]) * efunction[i]
            + c * (efunction[i - 1] + efunction[i + 1])
            + a * (efunction[i - 2] + efunction[i + 2]);
    }
    i = ngrid - 3;
    vtemp[i] = (b + vpot[i]) * efunction[i]
        + c * (efunction[i - 1] + efunction[i + 1])
        + a * efunction[i - 2];
    i = ngrid - 2;
    vtemp[i] = (b + vpot[i]) * efunction[i]
        + c * efunction[i + 1]
        + a * efunction[i - 2];
    //console.log(vtemp);
}

function invers(dx,vtemp,vpot,efunction,ngrid,nloop)
{
    let dtr = 1 / (dx * dx) / 2;
    let a = -dtr;
    let b = [];
    let bb = 2 * dtr;
    let c = -dtr;
    for (let k = 0; k < nloop; k++) {
        for (i = 0; i < ngrid; i++) {
            b[i] = (bb + vpot[i] - energy);
        }
        temp = b[0];
        efunction[0] = efunction[0] / temp;
        for (let j = 1; j < ngrid; j++) {
            temp1 = b[j];
            b[j] = c / temp;
            temp = temp1 - a * b[j];
            efunction[j] = (efunction[j] - a * efunction[j - 1]) / temp;
        }
        for (j = ngrid - 2; j > 1; j--) {
            efunction[j] = efunction[j] - b[j + 1] * efunction[j + 1];
        }
    }
}

function getEnergy()
{
    return initE;
}

// console.log(initE);

