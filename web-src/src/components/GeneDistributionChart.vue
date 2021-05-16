<style scoped>
</style>
<template>
<div id="dischart">
</div>
</template>
<script>
/* eslint-disable */
import {DistributionChart} from "@/amromicsjs";
// import SampleIGV from "@/components/Visualization/IGV";
export default {
  name: 'GeneDistributionChart',
  props: ['cluster_data'],
  data() {
    return {
      loading: false,

    };
  },
  mounted() {
    this.loading = true;
    var ctx=document.getElementById('dischart');
    //console.log(this.core_data);
    var dischart = new DistributionChart(ctx);
    dischart.load(this.cluster_data);

    dischart.draw();
    this.loading = false;
  },
  async created() {


  },
  method: {

  }
};
</script>
