library(shiny)

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(sf)
library(mapview)
library(leaflet)
library(htmltools)
library(BalancedSampling)

source("generate_Halton_pts.R")

reshapefile <- function(file_upload_data){
  ### rename uploaded files so it makes a sensible shapefile
  
  ## get the directory where shiny has dropped the uploads
  upload_dir = dirname(file_upload_data$datapath[1])
  message("Upload dir is ",upload_dir)
  ## need to rename to the original file name
  tos = file.path(upload_dir, file_upload_data$name)
  ## rename each of the uploaded files to the original file names
  file.rename(file_upload_data$datapath, tos)
  ## return the one with the .shp extention
  tos[grepl("\\.shp$", tos, ignore.case=TRUE)]
}

# read in study region shape file
#nc_sf <- st_read("data/mongolia.shp")

ui <- bootstrapPage(
  tags$style(type = "text/css", "html, body {width:100%;height:100%}"),
  leafletOutput("map", width = "100%", height = "100%"),
  absolutePanel(top = 10, right = 10,
                fileInput("dens_file",
                          "Locations of existing surveys",
                          multiple = TRUE),
                fileInput("survarea_file",
                          "Load new survey area",
                          multiple = TRUE),
                numericInput(inputId = "new_pts",
                             label = "Number of new survey points to generate",
                             min = 0,
                             value = 0),
                uiOutput("random_seed"),
                #numericInput(inputId = "seed",
                #             label = "Seed / survey ID",
                #             value = 1234),
                selectInput("basemap","Base map",
                            choices = list("Esri.WorldStreetMap", "Esri.WorldImagery", "OpenStreetMap"),
                            selected = "Esri.WorldStreetMap"),
                fluidRow(column(6, textInput(inputId = "strat_pr",
                                             label = "Sampling proportions per stratum")),
                         column(6, actionButton(inputId = "add_props", label = "Update"))),
                div(actionButton(inputId = "add_pts", label = "Add surveys"),
                    actionButton(inputId = "submit", label = "New seed"),
                    actionButton(inputId = "clear_pts", label = "Clear surveys"),
                    style = "margin:2px"),
                div(verbatimTextOutput("sb"), style = "margin:2px;width:200px;float:right")),
  absolutePanel(bottom = 10, left = 10,
                fluidRow(column(3, offset = 5, align="center",
                                downloadButton("downloadData", "Download survey points"))))
)

# Define server logic
server <- function(input, output, session) {
  
  nc_sf <- eventReactive(input$survarea_file, {
    shpfile = reshapefile(input$survarea_file)
    x <- st_read(shpfile)
    if(!("stratum_id" %in% names(x))){  x$stratum_id <- "all" }
    if(!is.factor(x$stratum_id)){ x$stratum_id <- factor(x$stratum_id) }
    x
  })
  
  # set colours for stratification variable
  pal <- reactive({colorFactor("YlOrRd", domain = nc_sf()$stratum_id)})
  
  micro_pts <- data.frame(lat = as.numeric(), lng = as.numeric())
  
  # stratum proportions
  strat_proportions <- eventReactive(c(nc_sf(), input$add_props), {
    # proportional to area covered
    sp <- nc_sf() %>% st_set_geometry(NULL) %>%
      count(stratum_id) %>%
      add_tally(n, name = "nn") %>%
      mutate(props = n / nn) %>%
      arrange(stratum_id) %>% select(props) %>% unlist()
    
    if(input$add_props > 0) sp <- eval(parse(text = input$strat_pr))
    
    # or specify these manually
    #sp <- c(0.1, 0.1, 0.8)
    
    sp
  })
  
  
  
  initial_seed <- 1234
  seed_value <- reactiveValues(seed = initial_seed)
  
  observeEvent(input$submit, {
    seed_value$seed <- round(runif(1, 1, 1e6))
  })
  
  # generate lots of Halton points
  
  # these points will be randomly generated across the study region; we need enough to make sure
  # there are *at least* the desired number of survey points in each stratum. There is no harm in
  # generating many more points that you need, because any contiguous subsequence is also spatially
  # well distributed, and the sequence is ordered (earlier elements surveyed first)
  
  pts_sf_ll <- eventReactive(c(input$seed, input$dens_file, input$survarea_file, input$add_props), {
    
    draw <- 10000
    set.seed(input$seed)
    my_seeds <- round(c(runif(1, 1, 1e6), runif(1, 1, 1e6)))
    pts <- generate_Halton_pts(n = draw, seeds = my_seeds, bases = c(2,3))
    
    # scale and shift Halton points to fit the bounding box of the study region
    bb <- st_bbox(nc_sf())
    scale <- c(bb$xmax, bb$ymax) - c(bb$xmin, bb$ymin)
    shift <- c(bb$xmin, bb$ymin)
    
    pts <- t(t(pts) * scale + shift)
    pts <- data.frame(pts)
    colnames(pts) <- c("X", "Y")
    
    # make the data frame an sf object
    pts_sf <- st_as_sf(pts, coords = c("X", "Y"), crs = st_crs(nc_sf())$proj4string)
    
    # remove any points that fall outside the study area
    pts_sf_in_studyarea <- lapply(st_intersects(pts_sf, nc_sf()), length) %>% unlist()
    pts_sf <- pts_sf[pts_sf_in_studyarea == 1, ]
    
    # remove any points that are too close to already surveyed points
    if (is.null(input$dens_file)) {
      pts_sf <- pts_sf
    } else {
      pts_sf_in_already_surveyed_area <- lapply(st_intersects(pts_sf, already_surveyed()$utm), length) %>% unlist()
      pts_sf <- pts_sf[pts_sf_in_already_surveyed_area == 0, ]
    }
    
    # work out what stratum each point falls into, and order the points sequentially both overall and within each stratum
    # (within-stratum order not necessary; global order used to determine survey ordering)
    pts_sf <- pts_sf %>% mutate(pt_id = 1:nrow(pts_sf),
                                stratum_id = nc_sf()$stratum_id[st_intersects(pts_sf, nc_sf()) %>% unlist()],
                                gridcell_id = nc_sf()$Id[st_intersects(pts_sf, nc_sf()) %>% unlist()]) %>%
      group_by(stratum_id) %>%
      mutate(ws_pt_id = row_number()) %>%
      ungroup() %>% cbind(st_coordinates(pts_sf))
    
    # add lat-long
    pts_sf_ll <- st_transform(pts_sf, crs = 4326)
    pts_sf_ll$lng <- st_coordinates(pts_sf_ll)[,1]
    pts_sf_ll$lat <- st_coordinates(pts_sf_ll)[,2]
    
    # add a variable indicating absolute order of adding points, this respects both Halton ordering
    # in pt_id and the desired number of points in each stratum
    obs_n_per_strat <- rep(0, length(strat_proportions()))
    strat_ids <- c()
    for(i in 1:nrow(pts_sf_ll)){
      exp_n_per_strat <- i * strat_proportions()
      this_strat <- which.min(obs_n_per_strat - exp_n_per_strat)
      strat_ids <- c(strat_ids, this_strat)
      obs_n_per_strat[this_strat] <- obs_n_per_strat[this_strat] + 1
    }
    order_df <- data.frame(ordered_id = 1:nrow(pts_sf_ll),
                           stratum_id = factor(levels(pts_sf_ll$stratum_id)[strat_ids],
                                               levels = levels(pts_sf_ll$stratum_id))) %>%
      group_by(stratum_id) %>%
      mutate(ws_pt_id = row_number()) %>% ungroup()
    
    
    pts_sf_ll <- left_join(pts_sf_ll, order_df, by = c("stratum_id", "ws_pt_id")) %>%
      arrange(ordered_id)
    
    # check on distributions across strata
    # table(pts_sf_ll$stratum_id[1:200])
    
    return(pts_sf_ll)
    
  })
  
  
  already_surveyed <- eventReactive(input$dens_file, {
    
    # x <- read_csv(input$dens_file$datapath)
    # 
    # # turn into sf object
    # already_sf <- st_as_sf(x, coords = c("lng", "lat"), crs = 4326)
    # # project for buffering
    # already_sf_utm <- st_transform(already_sf, crs = st_crs(nc_sf())$proj4string)
    # # add 10km buffer
    # already_sf_buffered <- st_buffer(already_sf_utm, dist = 10000)
    # # convert back to ll
    # already_sf_ll <- st_transform(already_sf_buffered, crs = 4326)
    
    shpfile = reshapefile(input$dens_file)
    already_sf_buffered <- st_read(shpfile)
    
    # convert back to ll
    already_sf_ll <- st_transform(already_sf_buffered, crs = 4326)
    
    return(list(latlong = already_sf_ll, utm = already_sf_buffered))
    
  })
  
  
  
  ######
  
  # want to keep track of (a) total number of points generated, (b) the last N points added,
  # to plot (a) as existing points and (b) as the latest set of points.
  
  # Reactive expression for the data subsetted to what the user selected
  existing_pts <- reactive({pts_sf_ll()[0,]})
  new_pts <- reactive({pts_sf_ll()[0,]})
  micro_pts <- reactive({pts_sf_ll()[0,]})
  
  # keep track of total points generated
  counter <- reactiveValues(countervalue = 0)
  
  observeEvent(input$add_pts, {
    counter$countervalue <- counter$countervalue + input$new_pts
  })
  
  # when counter changes, divide up pts_sf into (a) and (b) and the rest
  existing_pts <- eventReactive(c(input$add_pts, input$clear_pts), {
    
    n_existing <- max(0, counter$countervalue - input$new_pts)
    
    existing_pts <- pts_sf_ll() %>%
      filter(ordered_id <= n_existing)
    
    existing_pts <- st_as_sf(existing_pts, coords = c("lng", "lat"), crs = 4326, remove = FALSE)
    
    return(existing_pts)
    
  })
  
  new_pts <- eventReactive(existing_pts(), {
    
    n_new <- ifelse(counter$countervalue == 0, 0, input$new_pts)
    
    sampled_pts <- pts_sf_ll() %>%
      filter(ordered_id <= counter$countervalue)
    
    new_pts <- anti_join(sampled_pts, existing_pts() %>% st_set_geometry(NULL), by = "pt_id")
    
    return(new_pts)
    
  })
  
  # micro design points around an existing or new point
  micro_pts <- eventReactive(input$map_marker_click, {
    centroid_lat = input$map_marker_click$lat
    centroid_lng = input$map_marker_click$lng
    return(expand.grid(lat = seq(from = centroid_lat - 0.1, to = centroid_lat + 0.1, length.out = 5),
                       lng = seq(from = centroid_lng - 0.1, to = centroid_lng + 0.1, length.out = 5)))
  })
  
  # only include aspects of the map that won't need to change dynamically
  output$map <- renderLeaflet({
    leaflet() %>% addTiles() %>%
      flyToBounds(85, 30, 86, 53)
  })
  
  observe({
    leafletProxy("map") %>%
      clearMarkers() %>%
      clearGroup("maincells") %>%
      addPolygons(group = "maincells", data = st_transform(nc_sf(), crs = 4326), fillOpacity = 0.4, weight = 1, color = ~pal()(stratum_id)) %>%
      flyToBounds(min(pts_sf_ll()$lng), min(pts_sf_ll()$lat), max(pts_sf_ll()$lng), max(pts_sf_ll()$lat)) %>%
      addScaleBar(position = "bottomleft") %>%
      addLegend(layerId = "myleg", colors = pal()(levels(nc_sf()$stratum_id)), opacity = 0.7,
                position = "bottomright", labels = paste0(levels(nc_sf()$stratum_id), " (", 100*round(strat_proportions(), 2), "%)")) %>%
      # add inset map
      addMiniMap(
        tiles = providers$Esri.WorldStreetMap,
        position = 'bottomleft',
        width = 150, height = 150,
        zoomLevelOffset = -4,
        toggleDisplay = TRUE)
  })
  
  observe({
    leafletProxy("map") %>%
      clearMarkers() %>%
      addPolygons(data = already_surveyed()$latlong, fillOpacity = 0.2, weight = 0, color = "black")
  })
  
  observe({
    leafletProxy("map") %>%
      clearMarkers() %>%
      addCircleMarkers(data = new_pts(), layerId = ~pt_id, lng = ~lng, lat = ~lat, radius = 6, color = "red",
                       label = ~paste0("svy ", pt_id, " (lng ", round(lng,2), ", lat ", round(lat,2), ")")) %>%
      addCircleMarkers(data = existing_pts(), layerId = ~pt_id, radius = 3, lng = ~lng, lat = ~lat, color = "blue",
                       label = ~paste0("svy ", pt_id, " (lng ", round(lng,2), ", lat ", round(lat,2), ")")) %>%
      addCircleMarkers(data = st_as_sf(micro_pts(), coords = c("lng", "lat"), crs = 4326, remove = FALSE),
                       lng = ~lng, lat = ~lat, radius = 2, color = "yellow")
  })
  
  observe({
    leafletProxy("map") %>%
      clearTiles() %>%
      addTiles() %>%
      addProviderTiles(eval(parse(text=paste0("providers$",input$basemap))))
  })
  
  # measure of spatial balance
  output$sb <- renderText({
    N <- nrow(pts_sf_ll())
    n <- nrow(existing_pts()) + nrow(new_pts())
    p <- rep(n/N, N)
    X <- pts_sf_ll() %>% st_set_geometry(NULL) %>% select(X,Y) %>% as.matrix()
    s_exi <- existing_pts() %>% st_set_geometry(NULL) %>% select(pt_id) %>% unlist()
    s_new <- new_pts() %>% st_set_geometry(NULL) %>% select(pt_id) %>% unlist()
    s <- sort(c(s_exi, s_new))
    
    paste0("Spatial balance: ", round(BalancedSampling::sb(p, X, s), 3))
    
  })
  
  # remove all points on actionButton click
  observeEvent(c(input$clear_pts, input$seed, input$submit), {
    counter$countervalue <- 0
  })
  
  # OUTPUTS
  
  # generate the random seed input box
  output$random_seed <- renderUI({
    numericInput(inputId = "seed",
                 label = "Seed / survey ID",
                 value = seed_value$seed)
  })
  
  # downloadable csv of survey points
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("surveypts_seed", input$seed, ".csv")
    },
    content = function(file) {
      ex <- existing_pts() %>% st_set_geometry(NULL) %>% mutate(status = "existing") %>% select(pt_id, lat, lng, status, stratum_id, gridcell_id)
      new <- new_pts() %>% st_set_geometry(NULL) %>% mutate(status = "new") %>% select(pt_id, lat, lng, status, stratum_id, gridcell_id)
      file_to_write <- rbind(ex, new)
      write.csv(file_to_write, file, row.names = FALSE)
    }
  )
  
}

# Run the application
shinyApp(ui = ui, server = server)

