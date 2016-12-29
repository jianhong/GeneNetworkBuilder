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
    $.fn.cyExportbox = $.fn.cytoscapeExportbox = function(options){
      exportbox.apply(this, [options, $]);
      return this;
    };
    
    cytoscape('core', 'exportbox', function(options){
      var cy = this;
      exportbox.apply(cy.container(), [options, $]);
      return this;
    });
  };
  
  var defaults = {
    //parameters
    position     : 'bottom-right' //bottom-left, top-left, top-right
  };
  
  var exportbox = function(params, $){
    var options = $.extend(true, {}, defaults, params);
    var fn = params;
    var functions = {
      destroy: function(){
        var $this = $(this);
        var $sb = $this.find(".cy-exportbox");
        $sb.remove();
      },
      
      init: function(){
        return $(this).each(function(){
          var $container = $(this);
          
          var cy = $container.cytoscape('get');
          
          var $exportboxWrapper = $('<div class="cy-exportbox cy-exportbox-'+  
                                   options.position +'"></div>');
          $container.prepend($exportboxWrapper);
          var $exportbox = $('<div class="cy-exportbox-div"></div>');
          $exportboxWrapper.append($exportbox);
          
          $container.css('position', 'relative');
          $exportbox.css('position', 'absolute');
          
          var cyBtn = ['png', 'jpg', 'json'];
          
          function makeCyButton(type){
            var $btn = $('<button class="cy-exportbox-button" value=' +
                       type + '>' + type + '</button>');
            $btn.on('click', function(){
              var type=$(this).attr("value");
              var cy = $container.cytoscape('get');
              var a = document.createElement('a');
              if(type=='json'){
                a.href = 'data:application/json;charset=utf-8,'+ 
                          encodeURIComponent(JSON.stringify(cy.json()));
              }else{
                a.href = cy[type]({ full: true });
              }
              a.download = 'network.'+ type;
              function fireEvent(obj,evt){
                var fireOnThis = obj;
                var evObj;
                if( document.createEvent ) {
                  evObj = document.createEvent('MouseEvents');
                  evObj.initEvent( evt, true, false );
                  fireOnThis.dispatchEvent( evObj );
                } else if( document.createEventObject ) {
                  evObj = document.createEventObject();
                  fireOnThis.fireEvent( 'on' + evt, evObj );
                }
              }
              fireEvent(a, 'click');
            });
            $exportbox.append($btn);
          }
          
          cyBtn.forEach(makeCyButton);
        });
      }
    };
    if( functions[fn]){
      return functions[fn].apply(this, Array.prototype.slice.call(arguments, 1));
    } else if( typeof fn == 'object' || !fn ) {
      return functions.init.apply( this, arguments );
    } else {
      $.error("No such function `"+ fn +"` for jquery.cytoscapeExportbox");
    }
    return $(this);
  };
  
  if( typeof module !== 'undefined' && module.exports ){ // expose as a commonjs module
    module.exports = register;
  }

  if( typeof define !== 'undefined' && define.amd ){ // expose as an amd/requirejs module
    define('cytoscape-exportbox', function(){
      return register;
    });
  }

  if( typeof cytoscape !== 'undefined' && typeof jQuery !== 'undefined' ){ // expose to global cytoscape (i.e. window.cytoscape)
    register( cytoscape, jQuery );
  }

})();