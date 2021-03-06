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
        for(let j=i+1;j<n;j++){
            params[i][j]={'d':floyd[i][j],'k':K/floyd[i][j],'l':L*floyd[i][j]};
        }
    }
    return params;
}
function Kamada_Kawai(nodes, links){
    let n=nodes.length;
    let m=links.length;
    params = initgraph(nodes,links,n,m);
    for(let i=0;i<n;i++){
        nodes[i].x = Math.random() * 0.8 * width + 0.1 * width;
        nodes[i].y = Math.random() * 0.8 * height + 0.1 * height;
    }
    /*
    while(1){
        for(i=0;i<n;i++){
            
        }
    }
    */
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
    Kamada_Kawai(nodes, links);
    // 算法结束时间
    d2 = new Date()
    end = d2.getTime()

    // 保存图布局结果和花费时间为json格式，并按提交方式中重命名，提交时请注释这一部分代码
    var content = JSON.stringify({"time": end-begin, "nodes": nodes, "links": links});
    var blob = new Blob([content], { type: "text/plain;charset=utf-8" });
    //saveAs(blob, "save.json");
}

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

    let nodes_dict = {};
    for (i in nodes) {
        nodes_dict[nodes[i].id] = nodes[i]
    }

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
        .attr("fill", "steelblue")
        .on("mouseover", function (e, d) {// 鼠标移动到node上时显示text
            text
                .attr("display", function (f) {
                    if (f.id == d.id || f.weight > 40) {
                        return "null";
                    }
                    else {
                        return "none";
                    }
                })
        })
        .on("mouseout", function (e, d) {// 鼠标移出node后按条件判断是否显示text
            text
                .attr("display", function (f) {
                    if (f.weight > 40) {
                        return 'null';
                    }
                    else {
                        return 'none';
                    }
                })
        });

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

    // 图布局算法
    graph_layout_algorithm(nodes, links)

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

function main() {
    d3.json(data_file).then(function (DATA) {
        data = DATA;
        draw_graph();
    })
}

main()