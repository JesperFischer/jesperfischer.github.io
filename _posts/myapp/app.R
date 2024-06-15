
pacman::p_load(
  ordbetareg,
  dplyr,
  ggplot2,
  haven,
  brms,
  shinyWidgets,
  shiny,
  tidyr,
  stringr,
  extraDistr,
  ggh4x
)

#here::set_here()

source(here::here("_posts","myapp","utility_scripts.R"))

df = read.csv(here::here("_posts","myapp","raw_hrd.csv")) %>% filter(participant_id == "sub-0019" & Modality == "Extero")

dd = prep_data(df)

shiny::shinyApp(
  # Define UI
  ui <- navbarPage(
    "Parameters",
    tabPanel("JPDM", 
             sidebarLayout(
               sidebarPanel(
                 sliderInput("beta", "slope", min = -3, max = 3, value = 0.1, step = 0.1),
                 sliderInput("lapse", "lapse", min = -10, max = 0, value = -4, step = 0.01),
                 sliderInput("alpha", "alpha", min = -50, max = 50, value = 0, step = 0.1),
                 sliderInput("conf_int", "conf_int", min = -3, max = 3, value = -1, step = 0.1),
                 sliderInput("conf_beta", "conf_beta", min = -20, max = 5, value = -5, step = 0.1),
                 sliderInput("conf_shift", "conf_shift", min = -10, max = 10, value = 0, step = 0.1),
                 sliderInput("conf_cutzero", "conf_cutzero", min = -10, max = 10, value = 0, step = 0.1),
                 sliderInput("conf_cutone", "conf_cutone", min = -10, max = 10, value = 0, step = 0.1),
                 
                 sliderInput("rt_int", "rt_int", min = -3, max = 3, value = -0.5, step = 0.01),
                 sliderInput("rt_beta", "rt_beta", min = -3, max = 3, value = 1.5, step = 0.1),
                 sliderInput("rt_ndt", "rt_ndt", min = -10, max = 0, value = -3, step = 0.1),
                 sliderInput("rt_shift", "rt_shift", min = -10, max = 10, value = 0, step = 0.1)
                 
               ),
               mainPanel(
                 plotOutput("results", height = "900px")
               )
             )
    )
  ),
  
  # Define server
  server <- function(input, output) {
    
    
    
    output$results <- renderPlot({
      
      
      
      jpdm = function(x,beta,lapse,alpha,
                      conf_int,conf_beta,conf_prec,conf_shift,
                      rt_int,rt_beta,rt_sd,rt_ndt,rt_shift, 
                      id = 1){
        
        probab = function(x,lapse,beta,alpha){
          return((brms::inv_logit_scaled(lapse) / 2) + (1 - 2 *  (brms::inv_logit_scaled(lapse) / 2)) * (0.5+0.5 * pracma::erf((x-alpha) / (exp(beta) * sqrt(2)))))
          
        }
        prob = probab(x,lapse,beta,alpha)
        
        conf = brms::inv_logit_scaled(conf_int + conf_beta * (probab(x-conf_shift,lapse,beta,alpha) * (1-probab(x-conf_shift,lapse,beta,alpha))))
        
        rt = exp(rt_int + rt_beta * (probab(x-rt_shift,lapse,beta,alpha) * (1-probab(x-rt_shift,lapse,beta,alpha))))
        
        return(data.frame(prob = prob,rt = rt, conf  = conf,x = x, id = id))
      }
      # input = data.frame(beta = 1, lapse = 0.05, alpha = 5, conf_int = 1, conf_beta = -5, conf_prec = 10,conf_shift = -2,
      #                    rt_int = 2, rt_beta = 2, rt_sd = 1, rt_ndt = 1,rt_shift = 2, id = 1)  
      # 
      data <- jpdm(x = seq(-50,50,0.1), beta = input$beta, lapse = input$lapse, alpha = input$alpha,
                   conf_int = input$conf_int, conf_beta = input$conf_beta, conf_prec = input$conf_prec, input$conf_shift,
                   rt_int = input$rt_int, rt_beta = input$rt_beta, rt_sd = input$rt_sd, rt_ndt = input$rt_ndt, input$rt_shift, id = 1)
      
      
      #overlay data?
      
      dd = dd %>% mutate(conf = as.numeric(Confidence) / 100,
                         rt = as.numeric(rt), prob = y) %>% 
        pivot_longer(cols = c("prob","rt","conf"), names_to = "decision",values_to = "decision_value")
      
      
      data = data %>% pivot_longer(cols = c("prob","rt","conf"), names_to = "decision",values_to = "decision_value")
      
      plot1 = dd %>% 
        ggplot()+
        geom_point(data = dd, aes(x = x, y = decision_value))+
        geom_line(data = data,aes(x = x, y = decision_value))+
        facet_wrap(~decision, scales = "free",ncol = 1)+theme_classic()+
        theme(text = element_text(size = 24))
      
      
      
      plot1 = plot1+
        facetted_pos_scales(
          y = list(
            decision == "prob" ~ scale_y_continuous(limits = c(0, 1), breaks = c(0,0.2,0.4,0.6,0.8,1.0)),
            decision == "rt" ~ scale_y_continuous(limits = c(0, 8), breaks = c(0,2,4,6,8)),
            decision == "conf" ~ scale_y_continuous(limits = c(0, 1), breaks = c(0,0.2,0.4,0.6,0.8,1.0))
          )
        )
      
      
      plot1
    })
    
  },
  options = list(height = 1200, width = 1200)
  
)
