# kepler-mcmc.rb:
#   plots projections of the probability distributions for orbital
#   elements of G1 and G2.  probability distributions are determined
#   using the code in kepler-fit.py.
#
# Copyright 2015 Mike McCourt and Ann-Marie Madigan
#
# notes:
#   this uses the tioga plotting library for ruby (tioga.sf.net).
#   install it using "gem install tioga".
#
# usage:
#   create the plot in the appendix of McCourt & Madigan (2015) using
#   "tioga kepler-mcmc.rb -s 2"
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

    t.def_figure('g2-orbit-fit') do
      mnras_style
      enter_page_wide
      slurp_g2_table
      scatter_plots
    end

    t.def_figure('g1-orbit-fit') do
      mnras_style
      enter_page_wide
      slurp_g1_table
      scatter_plots
    end

    t.def_figure('g1-g2-orbit-fit') do
      mnras_style
      enter_page_special
      t.subplot('bottom_margin' => 2.5/9.0) do
        slurp_g1_table
        reverse_scatter_plots
      end
      t.subplot('top_margin' => 2.5/9.0) do
        slurp_g2_table
        scatter_plots
      end
    end

  end

  def enter_page
    mnras_style

    t.xlabel_shift = 2.0
    t.ylabel_shift = 1.75

    t.default_frame_left   = 0.18
    t.default_frame_right  = 0.96
    t.default_frame_top    = 0.96
    t.default_frame_bottom = 0.17

    t.default_page_width  = 72 * 2.3333

    t.default_page_height = t.default_page_width * \
      (t.default_frame_right - t.default_frame_left) / \
      (t.default_frame_top - t.default_frame_bottom)

    t.default_enter_page_function
  end

  def enter_page_wide
    mnras_style

    t.xlabel_shift = 2.0
    t.ylabel_shift = 1.75

    t.default_frame_left   = 0.06
    t.default_frame_right  = 0.98
    t.default_frame_top    = 0.97
    t.default_frame_bottom = 0.09

    t.default_page_width  = 72 * 7.0

    $nr = 7                     # number of plots in a row
    $margin = 0.0               # margin as a fraction of the page

    # size of a plot as a fraction of page width
    $plotsize = (1.0 - ($nr-1)*$margin)/$nr

    $nc = 7                     # number of plots in a column

    # total page height to width ratio
    $hfrac = ($nc * $plotsize + ($nc-1)*$margin)

    t.default_page_height = t.default_page_width * \
      (t.default_frame_right - t.default_frame_left) / \
      (t.default_frame_top - t.default_frame_bottom) * \
      $hfrac

    t.default_enter_page_function
  end

  def enter_page_special
    mnras_style

    t.xlabel_shift = 2.0
    t.ylabel_shift = 1.75

    t.default_frame_left   = 0.06
    t.default_frame_right  = 0.94
    t.default_frame_top    = 0.95
    t.default_frame_bottom = 0.05

    t.default_page_width  = 72 * 7.0

    $nr = 7                     # number of plots in a row
    $margin = 0.0               # margin as a fraction of the page

    # size of a plot as a fraction of page width
    $plotsize = (1.0 - ($nr-1)*$margin)/$nr

    $nc = 7                     # number of plots in a column

    # total page height to width ratio
    $hfrac = ($nc * $plotsize + ($nc-1)*$margin)

    t.default_page_height = t.default_page_width * \
      (t.default_frame_right - t.default_frame_left) / \
      (t.default_frame_top - t.default_frame_bottom) * \
      $hfrac * (9.0)/(9.0-2.5)

    t.default_enter_page_function
  end

  def slurp_g1_table
    ind, @a, @e, @i, @O, @o, @tp, @fp, @fv = \
    Dvector.fancy_read('chain-g1.dat') # trim

    @adata = [@a, 0.2,  0.9, '$a$ (arc-sec)']

    @edata = [@e, 0.8, 0.96, '$e$']

    @idata = [@i, 70.0, 77.0, '$i$']

    @Odata = [@O, 4.5, 35.0, '$\Omega$']

    @odata = [@o, 105.0, 130.0, '$\omega$']

    @tpdata = [@tp-2002, -1.5, 0.4, '$t_{\text{peri}}$ - 2002']

    @fpdata = [@fp, -5.0, -1.0, '$\ln f$']
  end

  def slurp_g2_table
    ind, @a, @e, @i, @O, @o, @tp, @fp, @fv = \
    Dvector.fancy_read('chain-g2.dat') # trim

    @adata = [@a, 0.5,  2.0, '$a$ (arc-sec)']

    @edata = [@e, 0.95, 1.0, '$e$']

    @idata = [@i, 56.0, 69.0, '$i$']

    @Odata = [@O, 0.0, 24.0, '$\Omega$']

    @odata = [@o, 87.0, 101.0, '$\omega$']

    @tpdata = [@tp-2014.25, -0.3, 0.2, '$t_{\text{peri}}$ - 2014.25']

    @fpdata = [@fp, -3.5, -2.0, '$\ln f$']
  end

  def make_histogram(xs, ys)
    xdata, xmin, xmax = xs
    ydata, ymin, ymax = ys

    n = 64
    tab = Dtable.new(n, n)

    xdata.each_index do |i|
      ix = (((xdata[i]-xmin) / (xmax-xmin)) * n).to_i
      iy = (((ydata[i]-ymin) / (ymax-ymin)) * n).to_i

      if (ix >= 0) and (ix < n) and (iy >= 0) and (iy < n)
        tab[iy,ix] += 1.0
      end
    end

    tab.div!(1.0 + tab.max).reverse_rows
  end

  def make_1D_histogram(xs)
    xdata, xmin, xmax = xs

    n = 64
    xvals = Dvector.new(n) {|i| xmin + (xmax-xmin)*i.to_f/(n-1)}
    yvals = Dvector.new(n)

    xdata.each do |x|
      ix = (((x-xmin) / (xmax-xmin)) * n).to_i
      if (ix >= 0) and (ix < n)
        yvals[ix] += 1
      end
    end

    [xvals, yvals.div!(yvals.max)]
  end

  def draw_1D_hist(xs, revearsed=false)
    xdata, xmin, xmax, xlabel = xs
    xvals, yvals = make_1D_histogram(xs)

    bins = Dvector[]
    tops = Dvector[]

    t.make_steps('xs' => xvals,
                 'ys' => yvals,
                 'dest_xs' => bins,
                 'dest_ys' => tops,
                 'x_first' => xmin,
                 'x_last'  => xmax,
                 'y_first' => 0.0,
                 'y_last'  => 0.0)


    t.ylabel_visible = false
    t.yaxis_type = AXIS_LINE_ONLY

    t.do_box_labels(nil, xlabel, nil)

    bounds = [xmin, xmax, 1.2, 0.0]
    if revearsed
      bounds = [xmin, xmax, 0.0, 1.2]
    end

    t.show_plot(bounds) do
      t.fill_color   = Green
      t.stroke_color = Black
      t.stroke_width = 0.5

      t.append_points_to_path(bins, tops)
      t.fill_and_stroke
    end
  end

  def draw_plot(xs, ys)
    xdata, xmin, xmax, xlabel = xs
    ydata, ymin, ymax, ylabel = ys

    t.do_box_labels(nil, xlabel, ylabel)

    hst = make_histogram([xdata, xmin, xmax], [ydata, ymin, ymax])

    t.show_plot([xmin, xmax, ymax, ymin]) do
      img = t.create_image_data(hst,
                                'min_value' => 0.0,
                                'max_value' => 1.0)

      t.show_image('ll' => [xmin, ymin],
                   'lr' => [xmax, ymin],
                   'ul' => [xmin, ymax],
                   'w'  => hst.num_cols,
                   'h'  => hst.num_rows,
                   'data' => img,
                   'color_space' => t.mellow_colormap)
    end

  end

  def scatter_plots
    space = ($nc.to_f) * ($plotsize + $margin)

    t.subplot('right_margin' => 6 * space/7,
              'left_margin'  => 0 * space/7) { a_column }

    t.ylabel_visible = false
    t.yaxis_type = AXIS_WITH_TICKS_ONLY

    t.subplot('right_margin' => 5 * space/7,
              'left_margin'  => 1 * space/7) { e_column }

    t.subplot('right_margin' => 4 * space/7,
              'left_margin'  => 2 * space/7) { i_column }

    t.subplot('right_margin' => 3 * space/7,
              'left_margin'  => 3 * space/7) { oo_column }

    t.subplot('right_margin' => 2 * space/7,
              'left_margin'  => 4 * space/7) { o_column }

    t.subplot('right_margin' => 1 * space/7,
              'left_margin'  => 5 * space/7) { tp_column }

    t.subplot('right_margin' => 0 * space/7,
              'left_margin'  => 6 * space/7) { fp_column }
  end

  def reverse_scatter_plots
    space = ($nc.to_f) * ($plotsize + $margin)

    t.ylabel_side = RIGHT
    t.yaxis_loc = RIGHT

    t.xlabel_side = TOP
    t.xaxis_loc = TOP

    t.subplot('right_margin' => 0 * space/7,
              'left_margin'  => 6 * space/7) { a_column(true) }

    t.ylabel_visible = false
    t.yaxis_type = AXIS_WITH_TICKS_ONLY

    t.subplot('right_margin' => 1 * space/7,
              'left_margin'  => 5 * space/7) { e_column(true) }

    t.subplot('right_margin' => 2 * space/7,
              'left_margin'  => 4 * space/7) { i_column(true) }

    t.subplot('right_margin' => 3 * space/7,
              'left_margin'  => 3 * space/7) { oo_column(true) }

    t.subplot('right_margin' => 4 * space/7,
              'left_margin'  => 2 * space/7) { o_column(true) }

    t.subplot('right_margin' => 5 * space/7,
              'left_margin'  => 1 * space/7) { tp_column(true) }

    t.subplot('right_margin' => 6 * space/7,
              'left_margin'  => 0 * space/7) { fp_column(true) }
  end

  def a_column(revearsed=false)
    space = ($nc.to_f) * ($plotsize + $margin)

    up   = 'top_margin'
    down = 'bottom_margin'

    if revearsed
      up   = 'bottom_margin'
      down = 'top_margin'
    end

    t.subplot(up   => 6 * space/7,
              down => 0 * space/7) { draw_plot(@adata, @fpdata) }

    t.xlabel_visible = false
    t.xaxis_type = AXIS_WITH_TICKS_ONLY

    t.subplot(up   => 5 * space/7,
              down => 1 * space/7) { draw_plot(@adata, @tpdata) }

    t.subplot(up   => 4 * space/7,
              down => 2 * space/7) { draw_plot(@adata, @odata) }

    t.subplot(up   => 3 * space/7,
              down => 3 * space/7) { draw_plot(@adata, @Odata) }

    t.subplot(up   => 2 * space/7,
              down => 4 * space/7) { draw_plot(@adata, @idata) }

    t.subplot(up   => 1 * space/7,
              down => 5 * space/7) { draw_plot(@adata, @edata) }

    t.subplot(up   => 0 * space/7,
              down => 6 * space/7) { draw_1D_hist(@adata, revearsed) }
  end

  def e_column(revearsed=false)
    space = ($nc.to_f) * ($plotsize + $margin)

    up   = 'top_margin'
    down = 'bottom_margin'

    if revearsed
      up   = 'bottom_margin'
      down = 'top_margin'
    end

    t.subplot(up   => 6 * space/7,
              down => 0 * space/7) { draw_plot(@edata, @fpdata) }

    t.xlabel_visible = false
    t.xaxis_type = AXIS_WITH_TICKS_ONLY

    t.subplot(up   => 5 * space/7,
              down => 1 * space/7) { draw_plot(@edata, @tpdata) }

    t.subplot(up   => 4 * space/7,
              down => 2 * space/7) { draw_plot(@edata, @odata) }

    t.subplot(up   => 3 * space/7,
              down => 3 * space/7) { draw_plot(@edata, @Odata) }

    t.subplot(up   => 2 * space/7,
              down => 4 * space/7) { draw_plot(@edata, @idata) }

    t.subplot(up   => 1 * space/7,
              down => 5 * space/7) { draw_1D_hist(@edata, revearsed) }
  end

  def i_column(revearsed=false)
    space = ($nc.to_f) * ($plotsize + $margin)

    up   = 'top_margin'
    down = 'bottom_margin'

    if revearsed
      up   = 'bottom_margin'
      down = 'top_margin'
    end

    t.subplot(up    => 6 * space/7,
              down => 0 * space/7) { draw_plot(@idata, @fpdata) }

    t.xlabel_visible = false
    t.xaxis_type = AXIS_WITH_TICKS_ONLY

    t.subplot(up    => 5 * space/7,
              down => 1 * space/7) { draw_plot(@idata, @tpdata) }

    t.subplot(up    => 4 * space/7,
              down => 2 * space/7) { draw_plot(@idata, @odata) }

    t.subplot(up    => 3 * space/7,
              down => 3 * space/7) { draw_plot(@idata, @Odata) }

    t.subplot(up    => 2 * space/7,
              down => 4 * space/7) { draw_1D_hist(@idata, revearsed) }
  end

  def oo_column(revearsed=false)
    space = ($nc.to_f) * ($plotsize + $margin)

    up   = 'top_margin'
    down = 'bottom_margin'

    if revearsed
      up   = 'bottom_margin'
      down = 'top_margin'
    end

    t.subplot(up    => 6 * space/7,
              down => 0 * space/7) { draw_plot(@Odata, @fpdata) }

    t.xlabel_visible = false
    t.xaxis_type = AXIS_WITH_TICKS_ONLY

    t.subplot(up    => 5 * space/7,
              down => 1 * space/7) { draw_plot(@Odata, @tpdata) }

    t.subplot(up    => 4 * space/7,
              down => 2 * space/7) { draw_plot(@Odata, @odata) }

    t.subplot(up    => 3 * space/7,
              down => 3 * space/7) { draw_1D_hist(@Odata, revearsed) }
  end

  def o_column(revearsed=false)
    space = ($nc.to_f) * ($plotsize + $margin)

    up   = 'top_margin'
    down = 'bottom_margin'

    if revearsed
      up   = 'bottom_margin'
      down = 'top_margin'
    end

    t.subplot(up    => 6 * space/7,
              down => 0 * space/7) { draw_plot(@odata, @fpdata) }

    t.xlabel_visible = false
    t.xaxis_type = AXIS_WITH_TICKS_ONLY

    t.subplot(up    => 5 * space/7,
              down => 1 * space/7) { draw_plot(@odata, @tpdata) }

    t.subplot(up    => 4 * space/7,
              down => 2 * space/7) { draw_1D_hist(@odata, revearsed) }
  end

  def tp_column(revearsed=false)
    space = ($nc.to_f) * ($plotsize + $margin)

    up   = 'top_margin'
    down = 'bottom_margin'

    if revearsed
      up   = 'bottom_margin'
      down = 'top_margin'
    end

    t.subplot(up    => 6 * space/7,
              down => 0 * space/7) { draw_plot(@tpdata, @fpdata) }

    t.xlabel_visible = false
    t.xaxis_type = AXIS_WITH_TICKS_ONLY

    t.subplot(up    => 5 * space/7,
              down => 1* space/7) { draw_1D_hist(@tpdata, revearsed) }
  end

  def fp_column(revearsed=false)
    space = ($nc.to_f) * ($plotsize + $margin)

    up   = 'top_margin'
    down = 'bottom_margin'

    if revearsed
      up   = 'bottom_margin'
      down = 'top_margin'
    end

    t.subplot(up    => 6 * space/7,
              down => 0 * space/7) { draw_1D_hist(@fpdata, revearsed) }
  end



end

MyPlots.new

# Local Variables:
#   compile-command: "tioga kepler-mcmc.rb -s"
#   coding: utf-8-unix
# End:
