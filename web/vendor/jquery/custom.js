// JavaScript Document

var elementThumbnailPhoto = 'ex-thumbnail-zoom-img';
var elementPhoto = 'ex-zoom-photo';

$('.' + elementThumbnailPhoto)
  // tile mouse actions
  .on('mouseover', function() {
    $(this).children('.' + elementPhoto).css({
      'transform': 'scale(' + $(this).attr('data-scale') + ')'
    });
  })
  .on('mouseout', function() {
    $(this).children('.' + elementPhoto).css({
      'transform': 'scale(1)'
    });
  })
  .on('mousemove', function(e) {
    $(this).children('.' + elementPhoto).css({
      'transform-origin': ((e.pageX - $(this).offset().left) / $(this).width()) * 100 + '% ' + ((e.pageY - $(this).offset().top) / $(this).height()) * 100 + '%'
    });
  })
  // tiles set up
  .each(function() {
    $(this)
      // add a photo container
      .append('<div class="' + elementPhoto + '"></div>')
      // set up a background image for each tile based on data-image attribute
      .children('.' + elementPhoto).css({
        'background-image': 'url(' + $(this).attr('data-image') + ')'
      });
  })