/*!
Copyright (c) Jianhong Ou
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the “Software”), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

;(function(){ 'use strict';
  //registers the extension on a cytoscape lib ref
  var register = function(cytoscape, $){
    if(!cytoscape){return;}
    $.fn.cySearchbox = $.fn.cytoscapeSearchbox = function(options){
      searchbox.apply(this, [options, $]);
      return this;
    };
    
    cytoscape('core', 'searchbox', function(options){
      var cy = this;
      searchbox.apply(cy.container(), [options, $]);
      return this;
    });
  };
  
  var defaults = {
    //parameters
    searchboxOnly  : false,
    blinkTime      : 500,
    backgroundColor: 'red',
    padding        : 30,
    searchIcon     : 'fa fa-search'
  };
  
  var searchbox = function(params, $){
    var options = $.extend(true, {}, defaults, params);
    var fn = params;
    var functions = {
      destroy: function(){
        var $this = $(this);
        var $sb = $this.find(".cy-searchbox");
        $sb.remove();
      },
      
      init: function(){
        return $(this).each(function(){
          var $container = $(this);
          
          var cy = $container.cytoscape('get');
          
          var $searchboxWrapper = $('<div class="cy-searchbox"></div>');
          $container.prepend($searchboxWrapper);
          var $searchbox = $('<div class="cy-searchbox-div"></div>');
          $searchboxWrapper.append($searchbox);
          
          $container.css('position', 'relative');
          $searchbox.css('position', 'absolute');
          
          if( options.searchboxOnly){
            $searchbox.addClass("cy-searchbox-searchbox-only");
          }
          
          var $querybox = $('<input type="text" class="cy-searchbox-querybox" autofocus placeholder="Search" />');
          $searchbox.append($querybox);
          
          var $submit = $('<div class="cy-searchbox-submit"><span class="icon ' + options.searchIcon + '"></span></div>');
          $searchbox.append($submit);
          
          function search(query){
            var nodes=[];
            var cy = $container.cytoscape("get");
            if(query!==""){
              cy.nodes().forEach(function(n){
                n.unselect();
              });
              nodes = cy.$(function(i, ele){
                if(ele.isNode()){
                  for(var key in ele.data()){
                    if(typeof ele.data(key) === "string"){
                      if(ele.data(key).toLowerCase().indexOf(query) > -1){
                        return true;
                      }
                    }
                  }
                }
              });
              if(nodes.length>0){
                nodes.forEach(function(n){
                  n.select();
                  for(var w=0; w<5; w++){
                    var bkcol = n.style('background-color');
                    n.animate({
                      style: { 'background-color': options.backgroundColor }
                      }, {
                        duration: options.blinkTime
                        }).delay( options.blinkTime ).animate({
                          style: { 'background-color': bkcol }
                          });
                          }
                });
                cy.fit(nodes, options.padding);
              }
            }
            return false;
          }
          
          $querybox.bind("keypress", function(e){
            var code = e.keyCode ? e.keyCode : e.which;
            if(code==13){
              search($(this).val().toLowerCase());
            }
          });
          
          $querybox.blur(function(){
            search($(this).val().toLowerCase());
          });
          
          $querybox.click(function(){
            $(this).select();
          });
          
          $submit.click(function(){
            search($querybox.val().toLowerCase());
          });
          
        });
      }
    };
    if( functions[fn]){
      return functions[fn].apply(this, Array.prototype.slice.call(arguments, 1));
    } else if( typeof fn == 'object' || !fn ) {
      return functions.init.apply( this, arguments );
    } else {
      $.error("No such function `"+ fn +"` for jquery.cytoscapeSearchbox");
    }
    return $(this);
  };
  
  if( typeof module !== 'undefined' && module.exports ){ // expose as a commonjs module
    module.exports = register;
  }

  if( typeof define !== 'undefined' && define.amd ){ // expose as an amd/requirejs module
    define('cytoscape-searchbox', function(){
      return register;
    });
  }

  if( typeof cytoscape !== 'undefined' && typeof jQuery !== 'undefined' ){ // expose to global cytoscape (i.e. window.cytoscape)
    register( cytoscape, jQuery );
  }

})();