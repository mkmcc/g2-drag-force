# parameter-space.rb:
#   makes figure 6 in McCourt & Madigan 2015.  computation is done in
#   the included python code drag-force-mla.py
#
# Copyright 2015 Mike McCourt and Ann-Marie Madigan
#
# notes:
#   this uses the tioga plotting library for ruby (tioga.sf.net).
#   install it using "gem install tioga".
#
# usage:
#   create the plot in the appendix of McCourt & Madigan (2015) using
#   "tioga parameter-space.rb -s 0"
#
require 'Tioga/FigureMaker'
require 'plot_styles.rb'
require 'Dobjects/Function'

class MyPlots

  include Math
  include Tioga
  include FigureConstants
  include MyPlotStyles

  def t
    @figure_maker
  end

  def initialize
    @figure_maker = FigureMaker.default

    t.tex_preview_preamble += "\n\t\\usepackage{mathtools}\n"
    t.tex_preview_preamble += "\n\t\\usepackage{textcomp}\n"
    t.tex_preview_preamble += "\n\t\\usepackage[charter]{mathdesign}\n"

    t.save_dir = 'plots'

    t.def_figure('scatter-plots-short') do
      mnras_style
      enter_page
      scatter_plots
    end

    slurp_table

    # tweak levels so that red = 1 sigma, gray = 2 sigma
    points = [0.0,   exp(-2.0), exp(-1.5), exp(-0.25), 1.0]
    colors = [White, LightGrey, DarkGrey,  FireBrick,  DarkRed]

    @my_colormap = \
    t.create_colormap('points' => points,
                      'Rs' => colors.map{|c| c[0]},
                      'Gs' => colors.map{|c| c[1]},
                      'Bs' => colors.map{|c| c[2]})
  end

  def enter_page
    mnras_style

    t.xlabel_shift = 2.0
    t.ylabel_shift = 1.75

    t.default_frame_left   = 0.06
    t.default_frame_right  = 0.98
    t.default_frame_top    = 0.96
    t.default_frame_bottom = 0.18

    t.default_page_width  = 72 * 7.0

    $nr = 3                     # number of plots in a row
    $margin = 0.1               # margin as a fraction of the page

    # size of a plot as a fraction of page width
    $plotsize = (1.0 - ($nr-1)*$margin)/$nr

    $nc = 1                     # number of plots in a column

    # total page height to width ratio
    $hfrac = ($nc * $plotsize + ($nc-1)*$margin)

    t.default_page_height = t.default_page_width * \
      (t.default_frame_right - t.default_frame_left) / \
      (t.default_frame_top - t.default_frame_bottom) * \
      $hfrac

    t.default_enter_page_function
  end

  def scatter_plots
    space = ($nr.to_f-1) * ($plotsize + $margin)

    t.subplot('right_margin' => space) { theta_phi_plot }

    t.subplot('left_margin'  => space/2,
              'right_margin' => space/2) { alpha_beta_plot }

    t.subplot('left_margin'  => space) { fkep_gf_plot }
  end

  def slurp_table
    puts "reading chain.dat..."

    @mc_ind,
    @mc_alpha,
    @mc_logbeta,
    @mc_fkep,
    @mc_theta,
    @mc_phi,
    @mc_loggf,
    @mc_raShift,
    @mc_decShift,
    @mc_vraShift,
    @mc_vdecShift,
    @mc_fp,
    @mc_fv = Dvector.fancy_read('chain-mla.dat')


    @mc_logbeta.div!(log(10.0))
    @mc_loggf.div!(log(10.0))
  end

  def make_histogram(xs, ys, title, n)
    # check if the data is cached
    cachefile = "plots/#{title}.hst"

    if File.exist?(cachefile)
      puts "  found cached data #{cachefile}... using that"
      #tab = Dtable.read(cachefile) # broken in new version of tioga?

      cols = Dvector.fancy_read(cachefile)
      m    = cols.length
      tab  = Dtable.new(m, m)
      cols.each_index {|i| tab.set_column(i, cols[i])}

      return tab
    end


    # read data if necessary and build the histogram
    if @mc_ind == nil
      slurp_table
    end

    puts "  building histogram from scratch"

    xdata, xmin, xmax = xs
    ydata, ymin, ymax = ys

    tab = Dtable.new(n, n)

    xdata.each_index do |i|
      ix = (((xdata[i]-xmin) / (xmax-xmin)) * n).to_i
      iy = (((ydata[i]-ymin) / (ymax-ymin)) * n).to_i

      if (ix >= 0) and (ix < n) and (iy >= 0) and (iy < n)
        tab[iy,ix] += 1.0
      end
    end

    tab = tab.div!(tab.max).reverse_rows


    # write the cache file so we don't need to remake this every time
    cols = Array.new

    (tab.num_cols).times {|i| cols << tab.column(i)}
    Dvector.write(cachefile, cols)

    return tab
  end

  def draw_2D_histogram(xs, ys, title, n=32)
    puts "drawing histogram #{title}..."

    xdata, xmin, xmax = xs
    ydata, ymin, ymax = ys

    hst = make_histogram([xdata, xmin, xmax], [ydata, ymin, ymax],
                         title, n)

    img = t.create_image_data(hst,
                              'min_value' => 0.0,
                              'max_value' => 1.0)

    t.show_image('ll' => [xmin, ymin],
                 'lr' => [xmax, ymin],
                 'ul' => [xmin, ymax],
                 'w'  => hst.num_cols,
                 'h'  => hst.num_rows,
                 'data' => img,
                 'color_space' => @my_colormap)

    dest_xs = Dvector.new; dest_ys = Dvector.new; gaps = Array.new
    xvals = Dvector.new(hst.num_cols) do |i|
      xmin + (xmax-xmin)*(0.5+i.to_f)/(hst.num_cols)
    end
    yvals = Dvector.new(hst.num_rows) do |i|
      ymin + (ymax-ymin)*(0.5+i.to_f)/(hst.num_rows)
    end
    dict = { 'dest_xs' => dest_xs,
             'dest_ys' => dest_ys,
             'gaps'    => gaps,
             'xs'      => xvals,
             'ys'      => yvals,
             'data'    => hst.reverse_rows }

    # plot 1- and 2-sigma contours
    #
    levels = [exp(-1.0), exp(-2.0)]
    levels.each do |level|
      dict['level'] = level
      t.make_contour(dict)
      t.append_points_with_gaps_to_path(dest_xs, dest_ys, gaps, false)
      t.stroke_width = 0.5
      t.stroke
    end
  end

  def theta_phi_plot
    t.do_box_labels(nil, '$\theta$', '$\phi$')

    t.xaxis_locations_for_major_ticks =\
    [0.0, PI/4, PI/2, 3*PI/4, PI]

    t.yaxis_locations_for_major_ticks =\
    [-PI, -PI/2, 0.0, PI/2, PI]

    t.xaxis_tick_labels =\
    ['0', '\pi/4', '\pi/2', '3\pi/4', '\pi']

    t.yaxis_tick_labels =\
    ['--\pi', '--\pi/2', '0', '\pi/2', '\pi']

    t.show_plot([0, PI, PI, -PI]) do
      draw_2D_histogram([@mc_theta, 0.0, PI], [@mc_phi, -PI, PI],
                        'theta-phi', 128)

      # draw a new x-axis
      degree = PI / 180
      t.xaxis_locations_for_major_ticks =\
      [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0].map{|x| x*degree}

      t.xaxis_tick_labels =\
      ['', '30$^{\circ}$', '60$^{\circ}$', '',
       '120$^{\circ}$', '150$^{\circ}$', '']

      t.show_axis({'from'          => [  0,0],
                   'to'            => [ PI,0],
                   'shift'         => -2.0,
                   'ticks_outside' => true,
                   'ticks_inside'  => true})


      # draw a new y-axis
      degree = PI / 180
      t.yaxis_locations_for_major_ticks =\
      [-180.0, -90.0, 0.0, 90.0, 180.0].map{|x| x*degree}

      t.yaxis_tick_labels =\
      ['', '--90$^{\circ}$', '', '90$^{\circ}$', '']

      t.show_axis({'from'          => [ PI/2, -PI],
                   'to'            => [ PI/2,  PI],
                   'shift'         => 1.25,
                   'ticks_outside' => true,
                   'ticks_inside'  => true})
    end
  end

  def fkep_gf_plot
    t.do_box_labels(nil,
                    '$f_{\text{kep}}$',
                    '$(L_{\text{cloud}} - R_{\text{cloud}})/R_{\text{cloud}}$')

    t.yaxis_log_values = true

    t.show_plot([0.0, 1.0, 2.0, 0.0]) do
      draw_2D_histogram([@mc_fkep, 0.0, 1.0],
                        [@mc_loggf, 0.0, 2.0],
                        'fkep-gf',
                        32)
    end
  end

  def alpha_beta_plot
    t.do_box_labels(nil, '$\alpha$', '$\beta$')

    t.yaxis_log_values = true

    t.show_plot([0.0, 1.2, 2.0, -2.0]) do
      draw_2D_histogram([@mc_alpha, 0.0, 1.2],
                        [@mc_logbeta, -2.0, 2.0],
                        'alpha-beta',
                        64)
    end
  end

end

MyPlots.new

# Local Variables:
#   compile-command: "tioga parameter-space.rb -s 0"
#   coding: utf-8-unix
# End:
