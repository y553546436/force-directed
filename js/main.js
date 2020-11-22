let _width = $(window).width();
let _height = $(window).height();
let width = _width;
let height = _height;

let data = null;
let data_file = './data/data.json';

function initgraph(nodes,links,n,m){
    name2id={};
    let inf=77777777;
    for(let i=0;i<n;i++){
        name2id[nodes[i].id] = i;
    }
    let floyd = []
    for(let i=0;i<n;i++){
        floyd[i]=[];
        for(j=0;j<n;j++)
            floyd[i][j]=inf;
    }
    for(let i=0;i<m;i++){
        links[i].from = name2id[links[i].source]
        links[i].to = name2id[links[i].target]
        floyd[links[i].from][links[i].to] = 1;
        floyd[links[i].to][links[i].from] = 1;
    }
    maxf = 0;
    for(let k=0;k<n;k++)
        for(let i=0;i<n;i++)
            for(let j=0;j<n;j++){
                if(i!=j&&j!=k&&i!=k)
                    floyd[i][j]=Math.min(floyd[i][j], floyd[i][k]+floyd[k][j]);
                if(k==n-1&&i<j){
                    maxf=Math.max(maxf,floyd[i][j]);
                }
            }
    L=height/maxf;
    K=233;
    params=[]
    for(let i=0;i<n;i++){
        params[i]=[]
        for(let j=0;j<n;j++){
            if(i!=j)
            params[i][j]={d:floyd[i][j],k:K/floyd[i][j],l:L*floyd[i][j]};
        }
    }
    /*
    for(let i=0;i<m;i++)
        if(links[i].from != links[i].to){
            params[links[i].from][links[i].to].k*=links[i].weight;
            params[links[i].to][links[i].from].k*=links[i].weight;
        }
        */
    return params;
}
function E_x(m, nodes, params) {
    let res = 0;
    let n=nodes.length;
    for(let i=0;i<n;i++)
        if(i!=m) {
            let k = nodes[m].x-nodes[i].x;
            let tmp=Math.sqrt(Math.pow(nodes[m].x-nodes[i].x,2)+Math.pow(nodes[m].y-nodes[i].y,2));
            k -= params[m][i].l*(nodes[m].x-nodes[i].x)/tmp;
            res += params[m][i].k*k;
        }
    return res;
}
function E_y(m, nodes, params) {
    let res = 0;
    let n=nodes.length;
    for(let i=0;i<n;i++)
        if(i!=m) {
            let k = nodes[m].y-nodes[i].y;
            let tmp=Math.sqrt(Math.pow(nodes[m].x-nodes[i].x,2)+Math.pow(nodes[m].y-nodes[i].y,2));
            k -= params[m][i].l*(nodes[m].y-nodes[i].y)/tmp;
            res += params[m][i].k*k;
        }
    return res;
}
function Delta(m, nodes, params) {
    return Math.sqrt(Math.pow(E_x(m,nodes,params),2)+Math.pow(E_y(m,nodes,params),2));
}
function E_xx(m,nodes,params) {
    let res = 0;
    let n=nodes.length;
    for(let i=0;i<n;i++)
        if(i!=m) {
            let k=1;
            let tmp=Math.pow(Math.pow(nodes[m].x-nodes[i].x,2)+Math.pow(nodes[m].y-nodes[i].y,2),1.5);
            k-=params[m][i].l*Math.pow(nodes[m].y-nodes[i].y,2)/tmp;
            res += params[m][i].k*k;
        }
    return res;
}
function E_xy(m,nodes,params) {
    let res = 0;
    let n=nodes.length;
    for(let i=0;i<n;i++)
        if(i!=m) {
            let k=0;
            let tmp=Math.pow(Math.pow(nodes[m].x-nodes[i].x,2)+Math.pow(nodes[m].y-nodes[i].y,2),1.5);
            k+=params[m][i].l*(nodes[m].y-nodes[i].y)*(nodes[m].x-nodes[i].x)/tmp;
            res += params[m][i].k*k;
        }
    return res;
}
function E_yx(m,nodes,params) {
    return E_xy(m,nodes,params);
}
function E_yy(m,nodes,params) {
    let res = 0;
    let n=nodes.length;
    for(let i=0;i<n;i++)
        if(i!=m) {
            let k=1;
            let tmp=Math.pow(Math.pow(nodes[m].x-nodes[i].x,2)+Math.pow(nodes[m].y-nodes[i].y,2),1.5);
            k-=params[m][i].l*Math.pow(nodes[m].x-nodes[i].x,2)/tmp;
            res += params[m][i].k*k;
        }
    return res;
}
function deltax(m,nodes,params) {
    return (E_y(m,nodes,params)*E_xy(m,nodes,params)-E_x(m,nodes,params)*E_yy(m,nodes,params))/(E_xx(m,nodes,params)*E_yy(m,nodes,params)-E_yx(m,nodes,params)*E_xy(m,nodes,params));
}
function deltay(m,nodes,params) {
    return (E_y(m,nodes,params)*E_xx(m,nodes,params)-E_x(m,nodes,params)*E_yx(m,nodes,params))/(E_xy(m,nodes,params)*E_yx(m,nodes,params)-E_yy(m,nodes,params)*E_xx(m,nodes,params));
}
function Kamada_Kawai(nodes, links){
    let n=nodes.length;
    let m=links.length;
    params = initgraph(nodes,links,n,m);
    for(let i=0;i<n;i++){
        nodes[i].x = Math.random() * 0.8 * width + 0.1 * width;
        nodes[i].y = Math.random() * 0.8 * height + 0.1 * height;
    }

    let loop_cnt=0;
    while(1){
        let max_Delta=0,k=0;
        for(i=0;i<n;i++){
            let Deltai=Delta(i,nodes,params);
            if(Deltai>max_Delta) {
                max_Delta=Deltai;
                k=i;
            }
        }
        let eps=1e5;
        if(max_Delta<eps) break;
        while(Delta(k,nodes,params)>eps) {
            delta_x=deltax(k,nodes,params);
            delta_y=deltay(k,nodes,params);
            nodes[k].x += delta_x;
            nodes[k].y += delta_y;
            loop_cnt++;
            if(loop_cnt==-1) {
                alert(k);
                alert(Delta(k,nodes,params));
                loop_cnt=0;
            }
        }
    }
}

let raw_l1=2,raw_l2=1e-6,raw_iters=500;
let lambda1=2,lambda2=1e-6,lambda3=0,iters=500;

function get_distance(a,b){
    return Math.sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}

function repulsive(nodes,i,j){
    let dist = get_distance(nodes[i],nodes[j]);
    let fs = -lambda1 * nodes[i].weight * nodes[j].weight/(dist);
    //fs = clip(fs);
    let dy = nodes[j].y - nodes[i].y, dx = nodes[j].x -nodes[i].x;
    return {x: fs*dx/dist, y: fs*dy/dist}
}

function attractive(nodes,i,j,w){
    let dist = get_distance(nodes[i],nodes[j]);
    //let len = lambda3 * nodes[i].weight * nodes[j].weight;
    //let fs = lambda2 * (dist*dist/len - len*len/dist);
    let fs = lambda2 * w * dist * dist;
    if(nodes[i].weight<=3 && nodes[j].weight>15) fs *= 10;
    if(nodes[j].weight<=3 && nodes[i].weight>15) fs *= 10;  
    //fs = clip(fs);
    let dy = nodes[j].y - nodes[i].y, dx = nodes[j].x -nodes[i].x;
    return {x: fs*dx/dist, y: fs*dy/dist}
}


function normalize(nodes) {
    n = nodes.length;
    let max_x=nodes[0].x,max_y=nodes[0].y,min_x=nodes[0].x,min_y=nodes[0].y;
    for(let i=1;i<n;++i) {
        max_x=Math.max(max_x,nodes[i].x);
        min_x=Math.min(min_x,nodes[i].x);
        max_y=Math.max(max_y,nodes[i].y);
        min_y=Math.min(min_y,nodes[i].y);
    }
    for(let i=0;i<n;++i) {
        nodes[i].x = (nodes[i].x-min_x)/(max_x-min_x)*0.7*width+0.05*width;
        nodes[i].y = (nodes[i].y-min_y)/(max_y-min_y)*0.8*height+0.1*height;
    }
}

/*
function boundforce(nodes,i) {
    let x=0,y=0;
    if(nodes[i].x>0.9*width) x=-(nodes[i].x-0.9*width);
    if(nodes[i].x<0.1*width) x=0.1*width-nodes[i].x;
    if(nodes[i].y>0.9*height) y=-(nodes[i].y-0.9*height);
    if(nodes[i].y<0.1*height) y=0.1*height-nodes[i].y;
    return {x : x*lambda3, y: y*lambda3};
}
*/
function clip(x){
    return x;
    x=Math.max(x, -1);
    x=Math.min(x, 1);
    return x;
}

function FR(nodes, links){
    let n=nodes.length;
    let m=links.length;
    name2id={};
    for(let i=0;i<n;i++){
        name2id[nodes[i].id] = i;
    }
    for(let i=0;i<m;i++){
        links[i].from = name2id[links[i].source]
        links[i].to = name2id[links[i].target]
    }
    for(let i=0;i<n;i++){
        nodes[i].x = Math.random() * 0.8 * width + 0.1 * width;
        nodes[i].y = Math.random() * 0.8 * height + 0.1 * height;
    }
    let tmp=1;
    for(let it=1;it<=iters;it++){
        //console.log(it)
        force=[]
        //tmp=0;
        for(let i=0;i<n;i++){
            force[i]={x:0,y:0};
            for(let j=0;j<n;j++){
                if(i!=j){
                    let fs = repulsive(nodes,i,j);
                    force[i].x+=fs.x;
                    force[i].y+=fs.y;
                    //console.log(fs.x, fs.y)
                    //tmp = Math.max(tmp, Math.max(Math.abs(fs.x),Math.abs(fs.y)))
                }
            }
        }
        //console.log('repulsive',tmp);
        //tmp=0;
        for(let e=0;e<m;e++){
            let x = links[e].from;
            let y = links[e].to;
            if(x==y) continue;
            let fs = attractive(nodes, x, y, links[e].weight);
            force[x].x+=fs.x;
            force[x].y+=fs.y;
            force[y].x-=fs.x;
            force[y].y-=fs.y;
            tmp = Math.max(tmp, Math.max(Math.abs(fs.x),Math.abs(fs.y)))
        }
        //console.log('attractive', tmp);
        /*
        for(let i=0;i<n;i++) {
            let fs = boundforce(nodes,i);
            force[i].x+=fs.x;
            force[i].y+=fs.y;
        }
        */
        cg =0;
        for(let i=0;i<n;i++){
            /*
            let fs = boundforce(nodes, i);
            force[i].x += fs.x;
            force[i].y += fs.y;
            */
            cg = Math.max(cg, Math.max(force[i].x, force[i].y));
            nodes[i].x += force[i].x;
            nodes[i].y += force[i].y;
        }
        //console.log(it, cg);
        //normalize(nodes);
        tmp*=0.9;
    }
    normalize(nodes);
    for(i in nodes){
        if(isNaN(nodes[i].x) || isNaN(nodes[i].y)){
            console.logs("Crash");
            redraw();
        }
    }
}
// 需要实现一个图布局算法，给出每个node的x,y属性
function graph_layout_algorithm(nodes, links) {
    // 算法开始时间
    d = new Date()
    begin = d.getTime()

    //这是一个随机布局，请在这一部分实现图布局算法
    /*
    for (i in nodes) {
        for (k = 0; k < 10000; k++) {
            nodes[i].x = Math.random() * 0.8 * width + 0.1 * width;
            nodes[i].y = Math.random() * 0.8 * height + 0.1 * height;
        }
    }
    */
    //Kamada_Kawai(nodes, links);
    FR(nodes,links);
    // 算法结束时间
    d2 = new Date()
    end = d2.getTime()

    // 保存图布局结果和花费时间为json格式，并按提交方式中重命名，提交时请注释这一部分代码
    var content = JSON.stringify({"time": end-begin, "nodes": nodes, "links": links});
    var blob = new Blob([content], { type: "text/plain;charset=utf-8" });
    //saveAs(blob, "save.json");
}
let nodes_dict = {};
let bgcolor = d3.color("rgba(240, 240, 240, 0.9)");

let col1 = d3.rgb(247,68,97);
let col3 = d3.rgb(87,96,105);
let col2 = d3.rgb(173,195,192);

function getc(w){
    if(w>40) return col1;
    else if(w>5) return col2;
    else return col3;
}
var drag = d3.drag()
        .on('drag', function (e, d) {
            d3.select(this).attr("cx", d.x = e.x ).attr("cy", d.y = e.y );
            d3.selectAll('line')
                .attr("x1", d => nodes_dict[d.source].x)
                .attr("y1", d => nodes_dict[d.source].y)
                .attr("x2", d => nodes_dict[d.target].x)
                .attr("y2", d => nodes_dict[d.target].y);
            d3.selectAll('text')
                .attr("x", d => d.x)
                .attr("y", d => d.y)
        });
function draw_graph() {
    let svg = d3.select('#container')
        .select('svg')
        .attr('width', width)
        .attr('height', height);

    // 数据格式
    // nodes = [{"id": 学校名称, "weight": 毕业学生数量}, ...]
    // links = [{"source": 毕业学校, "target": 任职学校, "weight": 人数}, ...]
    let links = data.links;
    let nodes = data.nodes;
    //len(nodes)=256 len(links)=846
    //console.log(links.length, nodes.length)
    // 图布局算法
    graph_layout_algorithm(nodes, links);

    let haslinks = []
    for(let i in nodes){
        haslinks[i]=[];
        for(let j in nodes)
            haslinks[i][j]=0;
    }
    for(let i in links){
        haslinks[links[i].from][links[i].to] = 1;
        haslinks[links[i].to][links[i].from] = 1;
    }
    for (i in nodes) {
        nodes_dict[nodes[i].id] = nodes[i]
    }

    for(i in nodes){
        nodes[i].rawcolor=getc(nodes[i].weight);
        let compute = d3.interpolate(bgcolor,nodes[i].rawcolor);
        nodes[i].fcolor=compute(0.2);
    }

    let clicking = false;

    // links
    let link = svg.append("g")
        .attr("stroke", "#999")
        .attr("stroke-opacity", 0.6)
        .selectAll("line")
        .data(links)
        .join("line")
        .attr("stroke-width", d => Math.sqrt(d.weight));

    // nodes
    let node = svg.append("g")
        .attr("stroke", "#fff")
        .attr("stroke-width", 0.5)
        .selectAll("circle")
        .data(nodes)
        .join("circle")
        .attr("r", d => Math.sqrt(d.weight) * 2 + 0.5)
        .attr("fill", d=>d.rawcolor)
        .on("mouseover", function (e, d) {// 鼠标移动到node上时显示text
            if(!clicking){
            text
                .attr("display", function (f) {
                    if (f.id == d.id || f.weight > 40) {
                        return "null";
                    }
                    else {
                        return "none";
                    }
                })
            }
        })
        .on("mouseout", function (e, d) {// 鼠标移出node后按条件判断是否显示text
            if(!clicking){
            text
                .attr("display", function (f) {
                    if (f.weight > 40) {
                        return 'null';
                    }
                    else {
                        return 'none';
                    }
                })
            }
        })
        .on("click", function (e, d){
            clicking = true;
            d3.selectAll("circle").attr("fill", (d2) =>{
                if(d==d2||haslinks[name2id[d.id]][name2id[d2.id]]) return d2.rawcolor;
                else return d2.fcolor;
            });
            d3.selectAll("line").style("visibility", (d2)=>{
                if(d2.source == d.id || d2.target == d.id) return "visible";
                else return "hidden";
            });
            text
                .attr("display", function(f){
                    if((f.id==d.id||haslinks[name2id[d.id]][name2id[f.id]])&&f.weight>5) return "null";
                    else return "none";
                })
        })
        .on("dblclick", function(e,d){
            clicking = false;
            d3.selectAll("circle").attr("fill", d=>d.rawcolor);
            d3.selectAll("line").style("visibility", "visible");
            text
                .attr("display", function (f) {
                    if (f.weight > 40) {
                        return 'null';
                    }
                    else {
                        return 'none';
                    }
                })
        })
        .call(drag);

    // 学校名称text，只显示满足条件的学校
    let text = svg.append("g")
        .selectAll("text")
        .data(nodes)
        .join("text")
        .text(d => d.id)
        .attr("display", function (d) {
            if (d.weight > 40) {
                return 'null';
            }
            else {
                return 'none';
            }
        });

    // 绘制links, nodes和text的位置
    link
        .attr("x1", d => nodes_dict[d.source].x)
        .attr("y1", d => nodes_dict[d.source].y)
        .attr("x2", d => nodes_dict[d.target].x)
        .attr("y2", d => nodes_dict[d.target].y);

    node
        .attr("cx", d => d.x)
        .attr("cy", d => d.y);
    text
        .attr("x", d => d.x)
        .attr("y", d => d.y)
    
}
function redraw(){
    d3.selectAll('svg > *').remove();
    draw_graph();
}

function slidechangerep(){
    v = document.getElementById('repr').value;
    lambda1 = raw_l1 * v;
    document.getElementById("bartext1").innerText="repulsive coefficient: " + v;

    d3.selectAll('svg > *').remove();
    draw_graph();
}

function slidechangeatt(){
    v = document.getElementById('attr').value;
    lambda2 = raw_l2 * v;
    document.getElementById("bartext2").innerText="attractive coefficient: " + v;

    d3.selectAll('svg > *').remove();
    draw_graph();
}

function slidechangeite(){
    v = document.getElementById('iter').value;
    iters = v;
    console.log(document.getElementById("bartext3").innerText)
    document.getElementById("bartext3").innerText="iterations: " + v;

    d3.selectAll('svg > *').remove();
    draw_graph();
}

function main() {
    d3.json(data_file).then(function (DATA) {
        data = DATA;
        draw_graph();
    })
}

main()
