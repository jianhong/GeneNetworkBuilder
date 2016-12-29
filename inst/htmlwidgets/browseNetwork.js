HTMLWidgets.widget({

  name: 'browseNetwork',

  type: 'output',

  factory: function(el, width, height) {
    var cy;
    return {

      renderValue: function(x) {

        //console.log(x);

        // TODO: code to render the widget, e.g.
        x.container = el;
        cy = cytoscape(x);
        //qtip
        cy.nodes().forEach(function(n){
          //console.log(n._private);
          var data = n._private.data;
          var g = data.label;
          var tipStr = "";
          // console.log(data);
          for(var key in data) {
            if(key === 'label'){
              tipStr += key + ': ' +
                '<a href="http://useast.ensembl.org/Multi/Search/Results?q=' +
                data[key] + '" target="_blank">' + data[key] + '</a><br/>';
            }else{
              tipStr += key + ": "+ data[key] + "<br/>";
            }
          }
          n.qtip({
            content: {text: tipStr,
                      title: g},
            position: {
              my: 'top center',
              at: 'bottom center'},
            style: {
              classes: 'qtip-light',
              tip: {
                width: 80,
                height: 8}
            }
          });
        });
        //panzoom
        cy.panzoom();
        //search box
        cy.searchbox();
        //export box
        cy.exportbox();
      },

      resize: function(width, height) {
        cy.resize();
      },

      cy: cy
    };
  }
});
