pipeline {
    agent any
    options {
        buildDiscarder(logRotator(numToKeepStr: '3'))
    }
    environment {
        IMAGE = 'mergeion2'
        NS = 'mergeion2'
        REGISTRY = '196229073436.dkr.ecr.eu-west-1.amazonaws.com'
        TAG = sh(returnStdout: true, script: "echo $BRANCH_NAME | sed -e 's/[A-Z]/\\L&/g' -e 's/[^a-z0-9._-]/./g'").trim()
        REGION = 'eu-west-1'
        NOT_CRAN = 'true'
    }
    stages {
        stage('Build Image') {
            agent {
                kubernetes {
                    yaml '''
                    apiVersion: v1
                    kind: Pod
                    spec:
                      containers:
                      - name: kaniko
                        image: gcr.io/kaniko-project/executor:v1.5.2-debug
                        env:
                        - name: AWS_SDK_LOAD_CONFIG
                          value: "true"
                        command:
                        - /kaniko/docker-credential-ecr-login
                        - get
                        tty: true
                        resources:
                          requests:
                              memory: "2048Mi"
                              ephemeral-storage: "2Gi"
                          limits:
                              memory: "4096Mi"
                              ephemeral-storage: "8Gi"
                        imagePullPolicy: Always
                      - name: aws-cli
                        image: amazon/aws-cli:2.1.33
                        command:
                        - cat
                        tty: true
                        resources:
                          requests:
                              memory: "100Mi"
                          limits:
                              memory: "2048Mi"'''
                    defaultContainer 'kaniko'
                }
            }
            steps {
                container('aws-cli') {
                    sh """aws --region ${env.REGION} ecr describe-repositories \
                    	--repository-names ${env.NS}/${env.IMAGE} \
                    	|| aws --region ${env.REGION} ecr create-repository \
                    	--repository-name ${env.NS}/${env.IMAGE}"""
                }
                container('kaniko') {
                    sh """/kaniko/executor \
                    	-v info \
                    	--context ${env.WORKSPACE} \
                    	--cache=true \
                    	--cache-ttl=8760h0m0s \
                    	--cache-repo ${env.REGISTRY}/${env.NS}/${env.IMAGE} \
                    	--cleanup \
                    	--destination ${env.REGISTRY}/${env.NS}/${env.IMAGE}:${env.TAG} \
                    	--registry-mirror ${env.REGISTRY}"""
                }
            }
        }
        stage('Packages') {
            agent {
                kubernetes {
                    yaml """
                    apiVersion: v1
                    kind: Pod
                    spec:
                      containers:
                      - name: r
                        image: ${env.REGISTRY}/${env.NS}/${env.IMAGE}:${env.TAG}
                        command: 
                        - cat
                        tty: true
                        imagePullPolicy: Always"""
                    defaultContainer 'r'
                }
            }
            stages {
                stage('MergeION') {
                    stages {
                        stage('Rcpp Compile Attributes') {
                            steps {
                                sh 'R -q -e \'if (requireNamespace("Rcpp", quietly = TRUE)) Rcpp::compileAttributes(".")\''
                            }
                        }
                        stage('Roxygen') {
                            steps {
                                sh 'R -q -e \'roxygen2::roxygenize(".")\''
                            }
                        }
                        stage('Build') {
                            steps {
                                sh 'R CMD build .'
                            }
                        }
                        stage('Check') {
                            steps {
                                script() {
                                    switch(sh(script: 'ls MergeION_*.tar.gz && R CMD check MergeION_*.tar.gz --no-manual', returnStatus: true)) {
                                        case 0: currentBuild.result = 'SUCCESS'
                                        default: currentBuild.result = 'FAILURE'; error('script exited with failure status')
                                    }
                                }
                            }
                        }
                        stage('Install') {
                            steps {
                                sh 'R -q -e \'install.packages(list.files(".", "MergeION_.*.tar.gz"), repos = NULL)\''
                            }
                        }
                        stage('Test and coverage') {
                            steps {
                                dir('.') {
                                    sh '''R -q -e \'code <- "testthat::test_package(\\"MergeION\\", reporter = testthat::MultiReporter$new(list(testthat::SummaryReporter$new(file = file.path(getwd(), \\"test-results.txt\\")), testthat::JunitReporter$new(file = file.path(getwd(), \\"results.xml\\")))))"\''''
                                   // packageCoverage <- covr::package_coverage(type = "none", code = code)
                                    // cat(readLines(file.path(getwd(), "test-results.txt")), sep = "\n")
                                    //covr::to_cobertura(packageCoverage)
                                }
                            }
                            /*post {
                                always {
                                    dir('.') {
                                        junit 'results.xml'
                                        cobertura autoUpdateHealth: false, autoUpdateStability: false, coberturaReportFile: 'cobertura.xml', conditionalCoverageTargets: '70, 0, 0', failUnhealthy: false, failUnstable: false, lineCoverageTargets: '80, 0, 0', maxNumberOfBuilds: 0, methodCoverageTargets: '80, 0, 0', onlyStable: false, sourceEncoding: 'ASCII', zoomCoverageChart: false
                                    }
                                }
                            } */
                        }
                    }
                }
                stage('Archive artifacts') {
                    steps {
                        archiveArtifacts artifacts: '*.tar.gz, *.pdf, **/00check.log', fingerprint: true
                    }
                }
            }
        }
    }
}

