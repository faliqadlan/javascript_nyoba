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

    // console.log(initE)

    
    //console.log(initE + 1e-16)
    singleEigen(xmin,dx,ngrid,v0,vpot,x,idx,efunction,iter,nloop,initE,leb,om, vtemp,potensial_function = harmonicPotential)
    //console.log(new Float32Array(initE))

    

    efunction.forEach(function(item) {
      var listItem = document.createElement('li');
      listItem.className = "nostyle";
      listItem.innerText = item;
      document.getElementById("ef").appendChild(listItem);
     });

}

// xmin = xmin.value;
// console.log(xmin,xmax,dx,initE,v0,omega);


    
function harmonicPotential(xmin,leb,om,x,v0 = 0)
{
    //$xp = [];
    xp = x - xmin - leb / 2;
    //console.log(xp, x, xmin, leb)
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
    //console.table(efunction)
    //Object.entries(efunction)
    vharm = vpot[ngrid - 1];
    // console.table(vharm)
    estart = initE;
    // console.table(estart)
    eps = 1e-5;
    if (estart > 0) eps = 1e-5;
    if (vharm > 1) eps = 1e-8;
    energy = initE;
    isig = 1;
    if (initE < 0) isig = -1;
    //console.table(isig)
    do {
        //console.log(initE);
        iter += iter;
        invers(dx,vtemp,vpot,efunction,ngrid,nloop);
        //console.log(efunction)
        sum = 0;
        for (let i = 0; i < ngrid; i++) {
            sum = sum + efunction[i] * efunction[i];
            //console.log(sum)
        }
        
        sum = Math.sqrt(sum * dx);
        //console.log(sum)
        for (i = 1; i < ngrid; i++) {
            efunction[i] = efunction[i] / sum;
            
        }
        //console.table(efunction)
        //console.log(x_tart)
        if (x_tart > 0) hamilton(dx,vtemp,vpot,efunction,ngrid);
        else hamilton5p(dx,vtemp,vpot,efunction,ngrid);
        //var_dump(global vtemp);
        //console.log(sum)
        sum = 0;
        //console.log(vtemp);
        for (let i = 1; i < ngrid - 1; i++) {
            sum = sum + vtemp[i] * efunction[i];
        }
        //console.log(sum)
        energy_new = sum * dx;
        //console.log(energy_new)
        delta = Math.abs((energy_new - energy) / energy);
        //console.log(eps,delta)
        energyN = energy_new;
        //console.log(iter)
        //console.log(energyN,initE)
        if (iter < 1) energyN = initE;
        energy = energyN;
        //console.log(energy,energyN)
    } while (delta <= eps);
    //console.log(initE, energy_new)
    initE = energy_new;

    // vtemp = [];
    let rsign = 1;
    if (x < 1e-8) rsign = (spaceship(efunction[1],0));
    //console.log(estart)
    ampl_max = 1;
    if (vharm <= 1) {
        if (estart > 0 && initE > 0) {
            ampl_max = 0;
            inode = 0;
            wb = efunction[ngrid - 1];
            //console.log(wb)
            i = ngrid;
            while (inode < 10) {
                i = i - 1;
                wf = efunction[i];
                //console.log(inode)
                if (i < 10) break;
                if (wf * wb < 0) inode = inode + 1;
                abswf = Math.abs(wf);
                //console.log(inode = inode + 1)
                if (abswf > ampl_max) ampl_max = abswf;
                wb = wf;
                //console.log(ampl_max)
            }
            //console.log(inode)
        }
    }
    //console.log(ampl_max)
    //console.log(efunction)
    ampl_max = ampl_max * rsign;
    for (let i = 0; i < ngrid; i++) {
        efunction[i] = efunction[i] / ampl_max;
    }
    //var_dump(efunction);
    //console.log(initE);
    // for (let i = 0; i < efunction.length; i++) {
    // return (
    //     cob.innerHTML = "hihia"
    // )
    

    // //     console.log(efunction[i]);
    // }
      //ef.innerHTML = efunction;

    //console.log(efunction)

    Enew.innerHTML = initE;

    // efunction.forEach(function(item) {
    //     var listItem = document.createElement('li');
    //     listItem.className = "nostyle";
    //     listItem.innerText = item;
    //     document.getElementById("ef").appendChild(listItem);
    //  });
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
    //console.log(vtemp);
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
    //console.table(efunction);
}

function invers(dx,vtemp,vpot,efunction,ngrid,nloop)
{
    //console.table(ngrid)
    let dtr = 1 / (dx * dx) / 2;
    let a = -dtr;
    let b = [];
    let bb = 2 * dtr;
    let c = -dtr;
    for (let k = 0; k < nloop; k++) {
        for (i = 0; i < ngrid; i++) {
            b[i] = (bb + vpot[i] - energy);
            // console.log(b[i])
        }
        //console.table(b)
        temp = b[0];
        //console.log(temp)
        efunction[0] = efunction[0] / temp;
        //console.table(efunction[0])
        //console.table(efunction[1])
        for (let j = 1; j < ngrid; j++) {
            temp1 = b[j];
            //console.table(temp1)
            b[j] = c / temp;
            // console.log(b[1])
            temp = temp1 - a * b[j];
            //console.log(temp)
            //console.table((efunction[j] - a * efunction[j - 1]) / temp)
            //console.log(efunction[1-1])
            efunction[j] = (efunction[j] - a * efunction[j - 1]) / temp;
            //console.log(efunction[1])
        }
        for (j = ngrid - 2; j > 1; j--) {
            
            efunction[j] = efunction[j] - b[j + 1] * efunction[j + 1];
            
        }
       // console.table(efunction)
    }
    // console.table(b)
}

function getEnergy()
{
    return initE;
}

// console.log(initE);

